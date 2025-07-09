// `timescale 1ns/1ps
/*  =======================  USER SWITCHES  ============================
 *  Edit if you want to change behaviour (no SV needed; plain defines)
 *  ------------------------------------------------------------------*/
`define STOP_WHEN_PERFECT   1   // 1 = break inner loop when correlation >= 0.95
`define DUMP_EVERY_BYTE     0   // 0 = only dump final snapshot (changed for performance)
`define DUMP_EVERY_ITER     0   // 0 = only dump final results (changed for performance)
`define CORRELATION_THRESH  0.95 // 0.95 correlation threshold
`define NUM_ITERATIONS      200  // Number of refinement iterations
/*  =================================================================== */

module tb_trng_attack_200_iter;

    /* ------------ DUT PARAMS (match RTL) ---------------------------- */
    parameter ADDR_WIDTH   = 14;   // 16 384 words → 128×128 image
    parameter DATA_WIDTH   = 52;
    parameter TRNG_A_WIDTH = 64;
    parameter TRNG_D_WIDTH = 32;

    /* ------------ expose user switches inside module ---------------- */
    parameter STOP_WHEN_PERFECT = `STOP_WHEN_PERFECT;
    parameter DUMP_EVERY_BYTE   = `DUMP_EVERY_BYTE;
    parameter DUMP_EVERY_ITER   = `DUMP_EVERY_ITER;
    parameter CORRELATION_THRESH = 950;  // Will be converted to real in task
    parameter NUM_ITERATIONS    = `NUM_ITERATIONS;

    /* ------------ clock / reset ------------------------------------ */
    reg clk;
    reg rst_n;
    initial begin
        clk = 0;
        rst_n = 0;
    end
    always #5 clk = ~clk;

    /* ------------ SRAM bus ----------------------------------------- */
    reg                    dcr, cs, we;
    reg  [ADDR_WIDTH-1:0]  addr;
    reg  [DATA_WIDTH-1:0]  wdata;
    reg  [TRNG_A_WIDTH-1:0] trng_a_in;
    reg  [TRNG_D_WIDTH-1:0] trng_d_in;
    wire [DATA_WIDTH-1:0]  rdata;
    wire                   ready;

    /* ------------ DUT instance ------------------------------------- */
    secure_sram_top #(
        .ADDR_WIDTH  (ADDR_WIDTH),
        .DATA_WIDTH  (DATA_WIDTH),
        .TRNG_A_WIDTH(TRNG_A_WIDTH),
        .TRNG_D_WIDTH(TRNG_D_WIDTH)
    ) dut (
        .clk       (clk),
        .rst_n     (rst_n),
        .dcr       (dcr),
        .cs        (cs),
        .we        (we),
        .addr      (addr),
        .wdata     (wdata),
        .trng_a_in (trng_a_in),
        .trng_d_in (trng_d_in),
        .rdata     (rdata),
        .ready     (ready)
    );

    /* ------------ reference image (8-bit pixels) ------------------- */
    reg [63:0] temp_pix [0:16383];   // Temporary for reading up to 44-bit hex values
    reg [7:0] plain_pix [0:16383];   // Final 8-bit pixels
    
    initial begin
        $readmemh("image_data.txt", temp_pix);
        // Convert to 8-bit pixels (take LSB)
        for (i = 0; i < 16384; i = i + 1) begin
            plain_pix[i] = temp_pix[i][7:0];
        end
        $display("Loaded image data: %0d pixels", 16384);
        $display("Sample pixels: %02X %02X %02X %02X %02X", 
                 plain_pix[0], plain_pix[1], plain_pix[2], plain_pix[3], plain_pix[4]);
    end

    /* ------------ secret & attacker keys --------------------------- */
    reg [63:0] true_trng_a;
    reg [31:0] true_trng_d;
    reg [31:0] bad_trng_d;
    initial begin
        true_trng_a = 64'hB4E7A1C9F82D5603;   // ground-truth address key
        true_trng_d = 32'h6F9C3E1A;           // ground-truth data key
        bad_trng_d  = 32'hDEADBEEF;           // attacker uses this
    end

    /* ------------ correlation calculation storage ------------------ */
    reg [7:0] decoded_pix [0:16383];   // store decoded pixels for correlation
    
    /* ------------ multi-iteration storage (EXPANDED) --------------- */
    reg [63:0] iteration_keys [0:199];   // store key from each iteration
    reg signed [31:0] iteration_scores [0:199];  // store correlation * 1000000 as integer
    
    /* ------------ locals ------------------------------------------- */
    integer i, val, best_val, iter_idx;
    reg signed [31:0] best_corr_int, cur_corr_int;  // correlation * 1000000
    real best_corr_real, cur_corr_real;             // for display only
    integer byte_idx, fh;
    reg [63:0] found_trng_a, cand;
    reg [255:0] fname;
    
    /* ------------ Progress tracking for long runs ------------------ */
    integer progress_milestone;
    real start_time, current_time;

    /* =================================================================
     *  MAIN FLOW
     * =================================================================*/
    initial begin
        $display("Starting 200-iteration TRNG attack simulation...");
        $display("Target key: %h", true_trng_a);
        $display("Progress will be reported every 10 iterations.");
        
        start_time = $realtime;
        
        #20 rst_n = 1;   // de-assert reset
        #50;  // Allow more settle time

        encode_reference_image();   // write obfuscated SRAM
        decode_reference_image();   // save full correct decode
        generate_comparison_case(); // save ideal comparison case

        $display("\n" + "="*60);
        $display("STARTING MULTI-ITERATION ATTACK (%0d iterations)", NUM_ITERATIONS);
        $display("="*60);

        /* ---------- Multi-iteration refinement loop ----------------- */
        found_trng_a = 64'h0;  // Start with all zeros

        for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
            // Progress reporting every 10 iterations
            if (iter_idx % 10 == 0) begin
                current_time = $realtime;
                $display("\n[PROGRESS] Iteration %0d/%0d (%.1f%%) - Time: %.2f ms", 
                         iter_idx, NUM_ITERATIONS, 
                         (iter_idx * 100.0) / NUM_ITERATIONS,
                         (current_time - start_time) / 1000000.0);
            end

            // Detailed progress for first few iterations only
            if (iter_idx < 3) begin
                $display("\n" + "="*40);
                $display("ITERATION %0d (starting key: %h)", iter_idx, found_trng_a);
                $display("="*40);
            end

            /* ---------- iterative brute-force per byte (MSB→LSB) ---- */
            for (byte_idx = 7; byte_idx >= 0; byte_idx = byte_idx - 1) begin
                if (iter_idx < 3) begin
                    $display("\n--- Iteration %0d, Searching byte %0d ---", iter_idx, byte_idx);
                end
                
                best_corr_int = -1000000;  // -1.0 * 1000000
                best_val = get_current_byte(found_trng_a, byte_idx);  // Start with current value

                for (val = 0; val < 256; val = val + 1) begin : SEARCH
                    cand = found_trng_a;
                    set_byte_value(cand, byte_idx, val);
                    pearson_score_candidate(cand, cur_corr_int);

                    if (cur_corr_int > best_corr_int) begin
                        best_corr_int = cur_corr_int;
                        best_val = val;
                    end
                    
                    // Early exit if correlation >= 0.95
                    if (STOP_WHEN_PERFECT && cur_corr_int >= 950000)
                        disable SEARCH;
                end

                set_byte_value(found_trng_a, byte_idx, best_val);
                best_corr_real = best_corr_int / 1000000.0;  // Convert back to real for display
                
                if (iter_idx < 3) begin
                    $display("iter%0d byte%0d chosen = 0x%02X  (correlation %0.6f)",
                             iter_idx, byte_idx, best_val[7:0], best_corr_real);
                end

                if (DUMP_EVERY_BYTE && iter_idx < 5) begin
                    dump_byte_snapshot(iter_idx, byte_idx);
                end
            end

            // Store iteration results
            iteration_keys[iter_idx] = found_trng_a;
            pearson_score_candidate(found_trng_a, iteration_scores[iter_idx]);
            
            if (iter_idx < 3) begin
                $display("\nIteration %0d complete:", iter_idx);
                $display("  Key: %h", found_trng_a);
                $display("  Score: %0.6f", iteration_scores[iter_idx] / 1000000.0);
            end

            if (DUMP_EVERY_ITER && iter_idx < 5) begin
                dump_iteration_snapshot(iter_idx);
            end
            
            // Check for early convergence
            if (found_trng_a == true_trng_a) begin
                $display("\n*** EXACT MATCH FOUND at iteration %0d! ***", iter_idx);
                $display("Key: %h", found_trng_a);
                $display("Score: %0.6f", iteration_scores[iter_idx] / 1000000.0);
                // Continue running to see if it stays stable
            end
            
            // Show high-scoring keys
            if (iteration_scores[iter_idx] > 800000) begin  // > 0.8 correlation
                $display("High score at iter %0d: %h (%.6f)", 
                         iter_idx, found_trng_a, iteration_scores[iter_idx] / 1000000.0);
            end
        end

        // Final summary
        current_time = $realtime;
        $display("\n" + "="*60);
        $display("MULTI-ITERATION ATTACK COMPLETED");
        $display("Total time: %.2f ms", (current_time - start_time) / 1000000.0);
        $display("="*60);
        
        analyze_convergence();
        dump_final_comparison();
        generate_summary_file();
        generate_convergence_analysis();
        
        $display("\nSimulation complete. Check output files for detailed results.");
        $finish;
    end

    /* =================================================================
     *  HELPER FUNCTIONS for byte manipulation (Pure Verilog)
     * =================================================================*/
    function [7:0] get_current_byte;
        input [63:0] key;
        input integer byte_idx;
        begin
            case (byte_idx)
                7: get_current_byte = key[63:56];
                6: get_current_byte = key[55:48];
                5: get_current_byte = key[47:40];
                4: get_current_byte = key[39:32];
                3: get_current_byte = key[31:24];
                2: get_current_byte = key[23:16];
                1: get_current_byte = key[15:8];
                0: get_current_byte = key[7:0];
                default: get_current_byte = 8'h00;
            endcase
        end
    endfunction

    task set_byte_value;
        inout [63:0] key;
        input integer byte_idx;
        input integer value;
        begin
            case (byte_idx)
                7: key[63:56] = value[7:0];
                6: key[55:48] = value[7:0];
                5: key[47:40] = value[7:0];
                4: key[39:32] = value[7:0];
                3: key[31:24] = value[7:0];
                2: key[23:16] = value[7:0];
                1: key[15:8]  = value[7:0];
                0: key[7:0]   = value[7:0];
            endcase
        end
    endtask

    /* =================================================================
     *  TASK: encode_reference_image
     *        write the plain pixels through DUT with TRUE keys
     * =================================================================*/
    task encode_reference_image;
        begin
            $display("Encoding reference image with true keys...");
            trng_a_in = true_trng_a;
            trng_d_in = true_trng_d;
            dcr = 1; 
            repeat (5) @(posedge clk); 
            dcr = 0; 
            repeat (5) @(posedge clk);

            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 1;
                addr  = i;
                wdata = {44'd0, plain_pix[i]};
                @(posedge clk);
                cs = 0; 
                we = 0;
                @(posedge clk);
            end
            repeat (10) @(posedge clk);
            $display("Reference image encoded successfully");
        end
    endtask

    /* =================================================================
     *  TASK: decode_reference_image
     *        true decode ─ used as correlation target
     * =================================================================*/
    task decode_reference_image;
        begin
            $display("Decoding reference image with true keys...");
            trng_a_in = true_trng_a;
            trng_d_in = true_trng_d;
            dcr = 1; 
            repeat (5) @(posedge clk); 
            dcr = 0; 
            repeat (5) @(posedge clk);

            fh = $fopen("reference_correct_decode.txt","w");
            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                repeat (2) @(posedge clk);
                $fwrite(fh, "%013X\n", rdata);
                cs = 0;
                @(posedge clk);
            end
            $fclose(fh);
            $display("Generated: reference_correct_decode.txt");
        end
    endtask

    /* =================================================================
     *  TASK: generate_comparison_case
     *        Generate ideal case: correct TRNG_A but wrong TRNG_D
     * =================================================================*/
    task generate_comparison_case;
        begin
            $display("Generating comparison case (correct TRNG_A, wrong TRNG_D)...");
            trng_a_in = true_trng_a;    // Use CORRECT address key
            trng_d_in = bad_trng_d;     // Use WRONG data key
            dcr = 1; 
            repeat (5) @(posedge clk); 
            dcr = 0; 
            repeat (5) @(posedge clk);

            fh = $fopen("correct_trnga_wrong_trngd.txt","w");
            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                repeat (2) @(posedge clk);
                $fwrite(fh, "%013X\n", rdata);
                cs = 0;
                @(posedge clk);
            end
            $fclose(fh);
            $display("Generated: correct_trnga_wrong_trngd.txt");
        end
    endtask

    /* =================================================================
     *  TASK: pearson_score_candidate
     *        Calculate Pearson correlation coefficient 
     *        Returns correlation * 1000000 as signed integer
     * =================================================================*/
    task pearson_score_candidate;
        input  [63:0] key;
        output reg signed [31:0] correlation_int;
        
        // Use 64-bit arithmetic to prevent overflow
        reg signed [63:0] sum_x, sum_y, sum_xy, sum_x2, sum_y2;
        reg signed [31:0] n;
        reg signed [63:0] numerator, denom_x, denom_y;
        real denom_real, correlation_real;
        
        // Debug variables
        integer debug_same_count, debug_zero_count, debug_max_decoded, debug_min_decoded;
        
        begin
            // Clear SRAM state and set new key
            trng_a_in = key;
            trng_d_in = bad_trng_d;
            dcr = 1; 
            repeat (5) @(posedge clk); 
            dcr = 0; 
            repeat (5) @(posedge clk);

            // Read all decoded pixels with proper timing
            debug_same_count = 0;
            debug_zero_count = 0;
            debug_max_decoded = 0;
            debug_min_decoded = 255;
            
            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                repeat (2) @(posedge clk);  // Extra cycles for stability
                decoded_pix[i] = rdata[7:0];  // Use LSB as decoded pixel
                cs = 0;
                @(posedge clk);
                
                // Debug statistics
                if (decoded_pix[i] == plain_pix[i]) debug_same_count = debug_same_count + 1;
                if (decoded_pix[i] == 0) debug_zero_count = debug_zero_count + 1;
                if (decoded_pix[i] > debug_max_decoded) debug_max_decoded = decoded_pix[i];
                if (decoded_pix[i] < debug_min_decoded) debug_min_decoded = decoded_pix[i];
            end
            
            // Calculate Pearson correlation coefficient using 64-bit arithmetic
            n = 16384;
            sum_x = 0; 
            sum_y = 0; 
            sum_xy = 0; 
            sum_x2 = 0; 
            sum_y2 = 0;
            
            // Calculate sums
            for (i = 0; i < n; i = i + 1) begin
                sum_x  = sum_x  + plain_pix[i];
                sum_y  = sum_y  + decoded_pix[i];
                sum_xy = sum_xy + (plain_pix[i] * decoded_pix[i]);
                sum_x2 = sum_x2 + (plain_pix[i] * plain_pix[i]);
                sum_y2 = sum_y2 + (decoded_pix[i] * decoded_pix[i]);
            end
            
            // Calculate correlation components
            numerator = (n * sum_xy) - (sum_x * sum_y);
            denom_x   = (n * sum_x2) - (sum_x * sum_x);
            denom_y   = (n * sum_y2) - (sum_y * sum_y);
            
            // Debug output for first iteration only
            if (iter_idx == 0 && byte_idx == 7 && val < 5) begin
                $display("DEBUG[key=%h]: same=%0d, zero=%0d, min=%0d, max=%0d", 
                         key, debug_same_count, debug_zero_count, debug_min_decoded, debug_max_decoded);
                $display("DEBUG: Plain[0:4]=%02X,%02X,%02X,%02X,%02X", 
                         plain_pix[0], plain_pix[1], plain_pix[2], plain_pix[3], plain_pix[4]);
                $display("DEBUG: Decoded[0:4]=%02X,%02X,%02X,%02X,%02X", 
                         decoded_pix[0], decoded_pix[1], decoded_pix[2], decoded_pix[3], decoded_pix[4]);
                $display("DEBUG: sum_x=%0d, sum_y=%0d, sum_xy=%0d", sum_x, sum_y, sum_xy);
                $display("DEBUG: denom_x=%0d, denom_y=%0d, numerator=%0d", denom_x, denom_y, numerator);
                
                // Show SRAM raw data for first few addresses
                $display("DEBUG: Raw SRAM data[0:4]=%h,%h,%h,%h,%h", 
                         rdata, rdata, rdata, rdata, rdata);
            end
            
            // Handle edge cases
            if (denom_x <= 0 || denom_y <= 0) begin
                correlation_int = -1000000;  // -1.0 * 1000000
                if (iter_idx == 0 && byte_idx == 7 && val < 5) begin
                    $display("DEBUG: Invalid correlation - negative variance");
                end
            end else begin
                // Calculate using real arithmetic for accuracy
                denom_real = $sqrt($itor(denom_x) * $itor(denom_y));
                correlation_real = $itor(numerator) / denom_real;
                
                // Clamp to valid range [-1.0, 1.0]
                if (correlation_real > 1.0) correlation_real = 1.0;
                if (correlation_real < -1.0) correlation_real = -1.0;
                
                // Convert to integer (* 1000000)
                correlation_int = $rtoi(correlation_real * 1000000.0);
                
                if (iter_idx == 0 && byte_idx == 7 && val < 5) begin
                    $display("DEBUG: Valid correlation = %.6f", correlation_real);
                end
            end
        end
    endtask

    /* =================================================================
     *  CONVERGENCE ANALYSIS
     * =================================================================*/
    task analyze_convergence;
        integer exact_matches, first_match_iter;
        reg signed [31:0] max_score;
        integer max_score_iter;
        integer stable_iterations;
        
        begin
            exact_matches = 0;
            first_match_iter = -1;
            max_score = -1000000;
            max_score_iter = -1;
            stable_iterations = 0;
            
            // Find statistics
            for (i = 0; i < NUM_ITERATIONS; i = i + 1) begin
                if (iteration_keys[i] == true_trng_a) begin
                    exact_matches = exact_matches + 1;
                    if (first_match_iter == -1) first_match_iter = i;
                end
                
                if (iteration_scores[i] > max_score) begin
                    max_score = iteration_scores[i];
                    max_score_iter = i;
                end
                
                // Check stability (key unchanged from previous iteration)
                if (i > 0 && iteration_keys[i] == iteration_keys[i-1]) begin
                    stable_iterations = stable_iterations + 1;
                end
            end
            
            $display("\nCONVERGENCE ANALYSIS:");
            $display("Ground-truth TRNG_A = %h", true_trng_a);
            $display("Exact matches: %0d/%0d iterations (%.1f%%)", 
                     exact_matches, NUM_ITERATIONS, 
                     (exact_matches * 100.0) / NUM_ITERATIONS);
            
            if (first_match_iter >= 0) begin
                $display("First exact match: iteration %0d", first_match_iter);
            end else begin
                $display("No exact matches found");
            end
            
            $display("Maximum score: %.6f at iteration %0d", max_score / 1000000.0, max_score_iter);
            $display("Best key: %h", iteration_keys[max_score_iter]);
            $display("Stable iterations: %0d/%0d (%.1f%%)", 
                     stable_iterations, NUM_ITERATIONS-1,
                     (stable_iterations * 100.0) / (NUM_ITERATIONS-1));
            
            // Show final 10 iterations
            $display("\nFinal 10 iterations:");
            for (i = NUM_ITERATIONS-10; i < NUM_ITERATIONS; i = i + 1) begin
                $display("  Iter %3d: %h (score: %.6f) %s", 
                         i, iteration_keys[i], iteration_scores[i] / 1000000.0,
                         (iteration_keys[i] == true_trng_a) ? "EXACT" : "");
            end
        end
    endtask

    /* =================================================================
     *  DUMP TASKS
     * =================================================================*/
    task dump_byte_snapshot;
        input integer iter_idx;
        input integer byte_idx;
        begin
            trng_a_in = found_trng_a;
            trng_d_in = bad_trng_d;
            dcr = 1; 
            repeat (5) @(posedge clk); 
            dcr = 0; 
            repeat (5) @(posedge clk);

            $sformat(fname, "iter%0d_byte%0d.txt", iter_idx, byte_idx);
            fh = $fopen(fname,"w");

            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                repeat (2) @(posedge clk);
                if (rdata === {DATA_WIDTH{1'bx}})
                    $fwrite(fh,"XXXXXXXXXXXXXXX\n");
                else
                    $fwrite(fh,"%013X\n", rdata);
                cs = 0;
                @(posedge clk);
            end
            $fclose(fh);
        end
    endtask

    task dump_iteration_snapshot;
        input integer iter_idx;
        begin
            trng_a_in = found_trng_a;
            trng_d_in = bad_trng_d;
            dcr = 1; 
            repeat (5) @(posedge clk); 
            dcr = 0; 
            repeat (5) @(posedge clk);

            $sformat(fname, "iteration_%0d_final.txt", iter_idx);
            fh = $fopen(fname,"w");

            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                repeat (2) @(posedge clk);
                if (rdata === {DATA_WIDTH{1'bx}})
                    $fwrite(fh,"XXXXXXXXXXXXXXX\n");
                else
                    $fwrite(fh,"%013X\n", rdata);
                cs = 0;
                @(posedge clk);
            end
            $fclose(fh);
            $display("   ↳ dumped %s", fname);
        end
    endtask

    task dump_final_comparison;
        integer sample_iterations[0:19]; // Sample every 10th iteration
        integer sample_idx;
        begin
            // Sample iterations: 0, 10, 20, ..., 190, 199
            for (i = 0; i < 19; i = i + 1) begin
                sample_iterations[i] = i * 10;
            end
            sample_iterations[19] = 199; // Include final iteration
            
            $display("Dumping sample iterations for analysis...");
            
            // Dump sampled iterations
            for (sample_idx = 0; sample_idx < 20; sample_idx = sample_idx + 1) begin
                iter_idx = sample_iterations[sample_idx];
                
                trng_a_in = iteration_keys[iter_idx];
                trng_d_in = bad_trng_d;
                dcr = 1; 
                repeat (5) @(posedge clk); 
                dcr = 0; 
                repeat (5) @(posedge clk);

                $sformat(fname, "final_iteration_%0d.txt", iter_idx);
                fh = $fopen(fname,"w");

                for (i = 0; i < 16384; i = i + 1) begin
                    @(posedge clk);
                    cs = 1; 
                    we = 0; 
                    addr = i;
                    @(posedge clk);
                    while (!ready) @(posedge clk);
                    repeat (2) @(posedge clk);
                    if (rdata === {DATA_WIDTH{1'bx}})
                        $fwrite(fh,"XXXXXXXXXXXXXXX\n");
                    else
                        $fwrite(fh,"%013X\n", rdata);
                    cs = 0;
                    @(posedge clk);
                end
                $fclose(fh);
            end
            $display("Generated sample iteration files for analysis");
        end
    endtask

    task generate_summary_file;
        begin
            fh = $fopen("attack_summary.txt","w");
            $fwrite(fh,"Multi-Iteration Attack Summary (200 iterations)\n");
            $fwrite(fh,"===============================================\n\n");
            $fwrite(fh,"Ground Truth TRNG_A: %h\n", true_trng_a);
            $fwrite(fh,"Number of iterations: %0d\n\n", NUM_ITERATIONS);
            
            // Write every 10th iteration + final few
            for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 10) begin
                $fwrite(fh,"Iteration %0d:\n", iter_idx);
                $fwrite(fh,"  Key:   %h\n", iteration_keys[iter_idx]);
                $fwrite(fh,"  Score: %0.6f\n", iteration_scores[iter_idx] / 1000000.0);
                $fwrite(fh,"  Match: %s\n\n", (iteration_keys[iter_idx] == true_trng_a) ? "EXACT" : "PARTIAL");
            end
            
            // Final iteration if not already included
            if ((NUM_ITERATIONS - 1) % 10 != 0) begin
                iter_idx = NUM_ITERATIONS - 1;
                $fwrite(fh,"Iteration %0d (FINAL):\n", iter_idx);
                $fwrite(fh,"  Key:   %h\n", iteration_keys[iter_idx]);
                $fwrite(fh,"  Score: %0.6f\n", iteration_scores[iter_idx] / 1000000.0);
                $fwrite(fh,"  Match: %s\n\n", (iteration_keys[iter_idx] == true_trng_a) ? "EXACT" : "PARTIAL");
            end
            
            $fclose(fh);
            $display("Generated: attack_summary.txt");
        end
    endtask

    task generate_convergence_analysis;
        begin
            fh = $fopen("convergence_analysis.txt","w");
            $fwrite(fh,"Convergence Analysis - 200 Iterations\n");
            $fwrite(fh,"====================================\n\n");
            
            // Write all iteration results for detailed analysis
            $fwrite(fh,"Iteration,Key,Score,Exact_Match\n");
            for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
                $fwrite(fh,"%0d,%h,%0.6f,%s\n", 
                        iter_idx, iteration_keys[iter_idx], iteration_scores[iter_idx] / 1000000.0,
                        (iteration_keys[iter_idx] == true_trng_a) ? "YES" : "NO");
            end
            
            $fclose(fh);
            $display("Generated: convergence_analysis.txt");
        end
    endtask

endmodule