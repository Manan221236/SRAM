// `timescale 1ns/1ps
/*  =======================  USER SWITCHES  ============================
 *  Edit if you want to change behaviour (no SV needed; plain defines)
 *  ------------------------------------------------------------------*/
`define STOP_WHEN_PERFECT   1   // 1 = break inner loop when correlation >= 0.95
`define DUMP_EVERY_BYTE     1   // 0 = only dump final snapshot
`define DUMP_EVERY_ITER     1   // 1 = dump after each iteration
`define CORRELATION_THRESH  0.95 // 0.95 correlation threshold
`define NUM_ITERATIONS      10  // Number of refinement iterations
/*  =================================================================== */

module tb_trngA_iterative;

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
    reg [7:0] plain_pix [0:16383];
    initial   $readmemh("image_data.txt", plain_pix);

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
    
    /* ------------ multi-iteration storage -------------------------- */
    reg [63:0] iteration_keys [0:9];   // store key from each iteration (10 iterations)
    real iteration_scores [0:9];       // store final score from each iteration (raw correlation)
    
    /* ------------ locals ------------------------------------------- */
    integer i, val, best_val, iter_idx;
    real best_corr, cur_corr;  // Use real for correlation values
    integer byte_idx, fh;
    reg [63:0] found_trng_a, cand;
    reg [255:0] fname;

    /* =================================================================
     *  MAIN FLOW
     * =================================================================*/
    initial begin
        #20 rst_n = 1;   // de-assert reset
        #20;

        encode_reference_image();   // write obfuscated SRAM
        decode_reference_image();   // save full correct decode
        generate_comparison_case(); // save ideal comparison case

        $display("\n" + "="*60);
        $display("STARTING MULTI-ITERATION ATTACK (%0d iterations)", NUM_ITERATIONS);
        $display("="*60);

        /* ---------- Multi-iteration refinement loop ----------------- */
        found_trng_a = 64'h0;  // Start with all zeros

        for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
            $display("\n" + "="*40);
            $display("ITERATION %0d (starting key: %h)", iter_idx, found_trng_a);
            $display("="*40);

            /* ---------- iterative brute-force per byte (MSB→LSB) ---- */
            for (byte_idx = 7; byte_idx >= 0; byte_idx = byte_idx - 1) begin
                $display("\n--- Iteration %0d, Searching byte %0d ---", iter_idx, byte_idx);
                best_corr = -1.0; 
                best_val = get_current_byte(found_trng_a, byte_idx);  // Start with current value

                for (val = 0; val < 256; val = val + 1) begin : SEARCH
                    cand = found_trng_a;
                    set_byte_value(cand, byte_idx, val);
                    pearson_score_candidate(cand, cur_corr);

                    if (cur_corr > best_corr) begin
                        best_corr = cur_corr;
                        best_val  = val;
                    end
                    if (STOP_WHEN_PERFECT && cur_corr >= 0.95)
                        disable SEARCH;               // early exit
                end

                set_byte_value(found_trng_a, byte_idx, best_val);
                $display("iter%0d byte%0d chosen = 0x%02X  (correlation %0.6f)",
                         iter_idx, byte_idx, best_val[7:0], best_corr);

                if (DUMP_EVERY_BYTE) dump_byte_snapshot(iter_idx, byte_idx);
            end

            // Store iteration results
            iteration_keys[iter_idx] = found_trng_a;
            pearson_score_candidate(found_trng_a, iteration_scores[iter_idx]);
            
            $display("\nIteration %0d complete:", iter_idx);
            $display("  Key: %h", found_trng_a);
            $display("  Score: %0.6f", iteration_scores[iter_idx]);

            if (DUMP_EVERY_ITER) dump_iteration_snapshot(iter_idx);
        end

        // Final summary
        $display("\n" + "="*60);
        $display("MULTI-ITERATION ATTACK SUMMARY");
        $display("="*60);
        $display("Ground-truth TRNG_A = %h", true_trng_a);
        for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
            $display("Iteration %0d result  = %h (score: %0.6f)", 
                     iter_idx, iteration_keys[iter_idx], iteration_scores[iter_idx]);
        end

        dump_final_comparison();
        generate_summary_file();
        $finish;
    end

    /* =================================================================
     *  HELPER FUNCTIONS for byte manipulation (Traditional Verilog)
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
            trng_a_in = true_trng_a;
            trng_d_in = true_trng_d;
            dcr = 1; 
            repeat (3) @(posedge clk); 
            dcr = 0; 
            @(posedge clk);

            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 1;
                addr  = i;
                wdata = {44'd0, plain_pix[i]};
                @(posedge clk);
                cs = 0; 
                we = 0;
            end
            repeat (5) @(posedge clk);
        end
    endtask

    /* =================================================================
     *  TASK: decode_reference_image
     *        true decode ─ used by Python as correlation target
     * =================================================================*/
    task decode_reference_image;
        begin
            trng_a_in = true_trng_a;
            trng_d_in = true_trng_d;
            dcr = 1; 
            repeat (3) @(posedge clk); 
            dcr = 0; 
            @(posedge clk);

            fh = $fopen("reference_correct_decode.txt","w");
            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                $fwrite(fh, "%013X\n", rdata);
                cs = 0;
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
            trng_a_in = true_trng_a;    // Use CORRECT address key
            trng_d_in = bad_trng_d;     // Use WRONG data key
            dcr = 1; 
            repeat (3) @(posedge clk); 
            dcr = 0; 
            @(posedge clk);

            fh = $fopen("correct_trnga_wrong_trngd.txt","w");
            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                $fwrite(fh, "%013X\n", rdata);
                cs = 0;
            end
            $fclose(fh);
            $display("Generated: correct_trnga_wrong_trngd.txt");
        end
    endtask

    /* =================================================================
     *  TASK: pearson_score_candidate
     *        Calculate Pearson correlation coefficient between 
     *        decoded pixels and original plaintext
     *        Returns correlation * 1000 as signed integer
     * =================================================================*/
    task pearson_score_candidate;
        input  [63:0] key;
        output real correlation;
        
        // Correlation calculation variables
        integer sum_x, sum_y, sum_xy, sum_x2, sum_y2;
        integer n;
        integer numerator, denom_x, denom_y;
        real denom_real;
        
        begin
            // First pass: read all decoded pixels
            trng_a_in = key;
            trng_d_in = bad_trng_d;
            dcr = 1; 
            repeat (3) @(posedge clk); 
            dcr = 0; 
            @(posedge clk);

            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                decoded_pix[i] = rdata[7:0];  // Use LSB as decoded pixel
                cs = 0;
            end
            
            // Calculate Pearson correlation coefficient
            n = 16384;
            sum_x = 0; 
            sum_y = 0; 
            sum_xy = 0; 
            sum_x2 = 0; 
            sum_y2 = 0;
            
            // First pass: calculate sums
            for (i = 0; i < n; i = i + 1) begin
                sum_x  = sum_x  + plain_pix[i];
                sum_y  = sum_y  + decoded_pix[i];
                sum_xy = sum_xy + (plain_pix[i] * decoded_pix[i]);
                sum_x2 = sum_x2 + (plain_pix[i] * plain_pix[i]);
                sum_y2 = sum_y2 + (decoded_pix[i] * decoded_pix[i]);
            end
            
            // Calculate correlation
            numerator = (n * sum_xy) - (sum_x * sum_y);
            denom_x   = (n * sum_x2) - (sum_x * sum_x);
            denom_y   = (n * sum_y2) - (sum_y * sum_y);
            
            // Handle edge cases
            if (denom_x <= 0 || denom_y <= 0) begin
                correlation = -1.0;  // Invalid correlation
            end else begin
                // Calculate denominator using real arithmetic
                denom_real = $sqrt($itor(denom_x) * $itor(denom_y));
                correlation = $itor(numerator) / denom_real;
                
                // Clamp to valid range [-1.0, 1.0]
                if (correlation > 1.0) correlation = 1.0;
                if (correlation < -1.0) correlation = -1.0;
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
            repeat (3) @(posedge clk); 
            dcr = 0; 
            @(posedge clk);

            $sformat(fname, "iter%0d_byte%0d.txt", iter_idx, byte_idx);
            fh = $fopen(fname,"w");

            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                if (rdata === {DATA_WIDTH{1'bx}})
                    $fwrite(fh,"XXXXXXXXXXXXXXX\n");
                else
                    $fwrite(fh,"%013X\n", rdata);
                cs = 0;
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
            repeat (3) @(posedge clk); 
            dcr = 0; 
            @(posedge clk);

            $sformat(fname, "iteration_%0d_final.txt", iter_idx);
            fh = $fopen(fname,"w");

            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1; 
                we = 0; 
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                if (rdata === {DATA_WIDTH{1'bx}})
                    $fwrite(fh,"XXXXXXXXXXXXXXX\n");
                else
                    $fwrite(fh,"%013X\n", rdata);
                cs = 0;
            end
            $fclose(fh);
            $display("   ↳ dumped %s", fname);
        end
    endtask

    task dump_final_comparison;
        begin
            // Dump each iteration's final result for Python analysis
            for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
                trng_a_in = iteration_keys[iter_idx];
                trng_d_in = bad_trng_d;
                dcr = 1; 
                repeat (3) @(posedge clk); 
                dcr = 0; 
                @(posedge clk);

                $sformat(fname, "final_iteration_%0d.txt", iter_idx);
                fh = $fopen(fname,"w");

                for (i = 0; i < 16384; i = i + 1) begin
                    @(posedge clk);
                    cs = 1; 
                    we = 0; 
                    addr = i;
                    @(posedge clk);
                    while (!ready) @(posedge clk);
                    @(posedge clk);
                    if (rdata === {DATA_WIDTH{1'bx}})
                        $fwrite(fh,"XXXXXXXXXXXXXXX\n");
                    else
                        $fwrite(fh,"%013X\n", rdata);
                    cs = 0;
                end
                $fclose(fh);
            end
            $display("Generated final iteration files for Python analysis");
        end
    endtask

    task generate_summary_file;
        begin
            fh = $fopen("attack_summary.txt","w");
            $fwrite(fh,"Multi-Iteration Attack Summary\n");
            $fwrite(fh,"==============================\n\n");
            $fwrite(fh,"Ground Truth TRNG_A: %h\n", true_trng_a);
            $fwrite(fh,"Number of iterations: %0d\n\n", NUM_ITERATIONS);
            
            for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
                $fwrite(fh,"Iteration %0d:\n", iter_idx);
                $fwrite(fh,"  Key:   %h\n", iteration_keys[iter_idx]);
                $fwrite(fh,"  Score: %0.6f\n", iteration_scores[iter_idx]);
                $fwrite(fh,"  Match: %s\n\n", (iteration_keys[iter_idx] == true_trng_a) ? "EXACT" : "PARTIAL");
            end
            $fclose(fh);
            $display("Generated: attack_summary.txt");
        end
    endtask

endmodule