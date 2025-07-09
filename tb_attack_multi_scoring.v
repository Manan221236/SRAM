`timescale 1ns / 1ps

module tb_attack_fixed;
    
    parameter NUM_ITERATIONS = 3;
    parameter NUM_TEST_PIXELS = 100;  // Test with smaller set first
    
    reg clk;
    reg rst_n;
    reg dcr;
    reg cs, we;
    reg [13:0] addr;
    reg [51:0] wdata;
    wire [51:0] rdata;
    wire ready;
    
    reg [63:0] trng_a_in;
    reg [31:0] trng_d_in;
    
    reg [7:0] plain_pix [0:NUM_TEST_PIXELS-1];
    reg [7:0] decoded_pix [0:NUM_TEST_PIXELS-1];
    reg [7:0] reference_pix [0:NUM_TEST_PIXELS-1];
    
    reg [63:0] true_trng_a;
    reg [31:0] true_trng_d;
    reg [63:0] best_trng_a;
    real best_correlation;
    
    integer i, byte_idx, val;
    real cur_correlation;
    integer best_val;
    
    // Module level variables
    integer current_byte_val;
    real effectiveness_ratio;
    real hamming_dist;
    
    secure_sram_top dut (
        .clk(clk),
        .rst_n(rst_n),
        .dcr(dcr),
        .trng_a_in(trng_a_in),
        .trng_d_in(trng_d_in),
        .cs(cs),
        .we(we),
        .addr(addr),
        .wdata(wdata),
        .rdata(rdata),
        .ready(ready)
    );
    
    initial clk = 0;
    always #5 clk = ~clk;
    
    initial begin
        $display("FIXED CORRELATION ATTACK");
        $display("========================");
        
        rst_n = 0;
        dcr = 0;
        cs = 0;
        we = 0;
        addr = 0;
        wdata = 0;
        trng_a_in = 0;
        trng_d_in = 0;
        best_correlation = -1.0;
        best_trng_a = 64'h0;
        
        true_trng_a = 64'hDEADBEEFCAFEBABE;
        true_trng_d = 32'h12345678;
        
        repeat (10) @(posedge clk);
        rst_n = 1;
        repeat (5) @(posedge clk);
        
        setup_attack;
        test_ground_truth;
        perform_attack;
        analyze_results;
        
        $finish;
    end

    task setup_attack;
        begin
            $display("\n=== SETUP ===");
            
            // Generate test pattern
            for (i = 0; i < NUM_TEST_PIXELS; i = i + 1) begin
                plain_pix[i] = (i * 7 + 13) % 256;  // Structured pattern
            end
            
            $display("Generated %0d test pixels", NUM_TEST_PIXELS);
            $display("Sample: %0d %0d %0d %0d %0d", 
                     plain_pix[0], plain_pix[1], plain_pix[2], plain_pix[3], plain_pix[4]);
            
            // Encrypt with true keys
            encrypt_data;
            
            // Create reference by reading with true keys
            create_reference;
        end
    endtask
    
    task encrypt_data;
        begin
            $display("Encrypting data with true keys...");
            
            complete_system_reset;
            
            trng_a_in = true_trng_a;
            trng_d_in = true_trng_d;
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (5) @(posedge clk);

            // Write data (gets encrypted)
            for (i = 0; i < NUM_TEST_PIXELS; i = i + 1) begin
                @(posedge clk);
                cs = 1;
                we = 1;
                addr = i;
                wdata = {43'd0, plain_pix[i]};
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                cs = 0;
                we = 0;
                
                // Small delay between operations
                repeat (2) @(posedge clk);
            end
            
            $display("Data encrypted and stored");
            complete_system_reset;
        end
    endtask
    
    task create_reference;
        integer matches;
        begin
            $display("Creating reference by reading with true keys...");
            
            complete_system_reset;
            
            trng_a_in = true_trng_a;
            trng_d_in = true_trng_d;
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (5) @(posedge clk);

            // Read back with true keys (should give original)
            for (i = 0; i < NUM_TEST_PIXELS; i = i + 1) begin
                @(posedge clk);
                cs = 1;
                we = 0;
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                reference_pix[i] = rdata[7:0];
                cs = 0;
                
                repeat (2) @(posedge clk);
            end
            
            $display("Reference created. Sample: %0d %0d %0d %0d %0d", 
                     reference_pix[0], reference_pix[1], reference_pix[2], reference_pix[3], reference_pix[4]);
            
            // Check how many match
            matches = 0;
            for (i = 0; i < NUM_TEST_PIXELS; i = i + 1) begin
                if (reference_pix[i] == plain_pix[i]) matches = matches + 1;
            end
            
            $display("Reference check: %0d/%0d pixels match original (%0.1f%%)", 
                     matches, NUM_TEST_PIXELS, (matches * 100.0) / NUM_TEST_PIXELS);
            
            if (matches > NUM_TEST_PIXELS * 0.9) begin
                $display("GOOD: Reference matches original - system working");
            end else begin
                $display("ISSUE: Reference doesn't match - but continuing attack");
            end
            
            complete_system_reset;
        end
    endtask

    task test_ground_truth;
        begin
            $display("\n=== GROUND TRUTH TEST ===");
            
            // Test true key
            calculate_correlation(true_trng_a, cur_correlation);
            $display("True TRNG_A correlation: %0.6f", cur_correlation);
            
            // Test wrong keys
            calculate_correlation(64'h0, cur_correlation);
            $display("Zero TRNG_A correlation: %0.6f", cur_correlation);
            
            calculate_correlation(64'hFFFFFFFFFFFFFFFF, cur_correlation);
            $display("Ones TRNG_A correlation: %0.6f", cur_correlation);
            
            if (cur_correlation > 0.5) begin
                $display("WARNING: Random key gives high correlation - attack may not work");
            end
        end
    endtask

    task perform_attack;
        begin
            $display("\n=== PERFORMING ATTACK ===");
            $display("Starting from all zeros and optimizing byte-by-byte...");
            
            best_trng_a = 64'h0000000000000000;
            calculate_correlation(best_trng_a, best_correlation);
            $display("Initial key: %h, correlation: %0.6f", best_trng_a, best_correlation);
            
            // Optimize each byte
            for (byte_idx = 7; byte_idx >= 0; byte_idx = byte_idx - 1) begin
                optimize_byte(byte_idx);
            end
            
            $display("\nFinal key: %h", best_trng_a);
            $display("Final correlation: %0.6f", best_correlation);
        end
    endtask
    
    task optimize_byte;
        input integer byte_pos;
        
        real best_byte_corr, test_corr;
        integer best_byte_val;
        reg [63:0] test_key;
        
        begin
            $display("\n--- Optimizing Byte %0d ---", byte_pos);
            
            best_byte_corr = best_correlation;
            current_byte_val = get_current_byte(best_trng_a, byte_pos);
            best_byte_val = current_byte_val;
            
            $display("Current value: 0x%02X, testing all 256 values...", current_byte_val);
            
            // Try all 256 values
            for (val = 0; val < 256; val = val + 1) begin
                test_key = best_trng_a;
                set_byte_value(test_key, byte_pos, val);
                
                calculate_correlation(test_key, test_corr);
                
                if (test_corr > best_byte_corr) begin
                    best_byte_corr = test_corr;
                    best_byte_val = val;
                    $display("  New best: 0x%02X -> correlation %0.6f", val, test_corr);
                end
                
                // Progress indicator
                if (val % 64 == 63) begin
                    $display("  Progress: %0d/256 tested, best so far: 0x%02X (%0.6f)", 
                             val + 1, best_byte_val, best_byte_corr);
                end
            end
            
            // Apply best value
            set_byte_value(best_trng_a, byte_pos, best_byte_val);
            best_correlation = best_byte_corr;
            
            $display("BYTE %0d RESULT: 0x%02X (correlation: %0.6f)", 
                     byte_pos, best_byte_val, best_byte_corr);
            $display("Current key: %h", best_trng_a);
        end
    endtask

    task calculate_correlation;
        input [63:0] test_key;
        output real correlation;
        
        integer sum_x, sum_y, sum_xy, sum_x2, sum_y2;
        integer n, numerator, denom_x, denom_y;
        real denom_real;
        integer matches;
        
        begin
            complete_system_reset;
            
            trng_a_in = test_key;
            trng_d_in = true_trng_d;
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (5) @(posedge clk);

            // Read decoded pixels
            for (i = 0; i < NUM_TEST_PIXELS; i = i + 1) begin
                @(posedge clk);
                cs = 1;
                we = 0;
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                decoded_pix[i] = rdata[7:0];
                cs = 0;
                
                repeat (2) @(posedge clk);
            end

            // Calculate Pearson correlation
            sum_x = 0; sum_y = 0; sum_xy = 0; sum_x2 = 0; sum_y2 = 0;
            n = NUM_TEST_PIXELS;
            
            for (i = 0; i < NUM_TEST_PIXELS; i = i + 1) begin
                sum_x = sum_x + reference_pix[i];
                sum_y = sum_y + decoded_pix[i];
                sum_xy = sum_xy + (reference_pix[i] * decoded_pix[i]);
                sum_x2 = sum_x2 + (reference_pix[i] * reference_pix[i]);
                sum_y2 = sum_y2 + (decoded_pix[i] * decoded_pix[i]);
            end
            
            numerator = (n * sum_xy) - (sum_x * sum_y);
            denom_x = (n * sum_x2) - (sum_x * sum_x);
            denom_y = (n * sum_y2) - (sum_y * sum_y);
            
            if (denom_x <= 0 || denom_y <= 0) begin
                correlation = 0.0;  // Avoid NaN
            end else begin
                denom_real = $sqrt($itor(denom_x) * $itor(denom_y));
                correlation = $itor(numerator) / denom_real;
                
                if (correlation > 1.0) correlation = 1.0;
                if (correlation < -1.0) correlation = -1.0;
            end
            
            // Count exact matches for debugging
            matches = 0;
            for (i = 0; i < NUM_TEST_PIXELS; i = i + 1) begin
                if (decoded_pix[i] == reference_pix[i]) matches = matches + 1;
            end
            
            // Debug key cases
            if (test_key == true_trng_a || test_key == 64'h0 || correlation > 0.8) begin
                $display("    Key %h: %0d/%0d exact matches, correlation=%0.6f", 
                         test_key, matches, NUM_TEST_PIXELS, correlation);
            end
        end
    endtask

    task analyze_results;
        real hamming_distance;
        integer correct_bytes;
        begin
            $display("\n=== RESULTS ===");
            $display("Best found TRNG_A: %h", best_trng_a);
            $display("True TRNG_A:       %h", true_trng_a);
            $display("Best correlation:  %0.6f", best_correlation);
            
            hamming_distance = calculate_hamming_distance(best_trng_a, true_trng_a);
            correct_bytes = count_correct_bytes(best_trng_a, true_trng_a);
            
            $display("Hamming distance:  %0.1f bits", hamming_distance);
            $display("Correct bytes:     %0d/8", correct_bytes);
            
            if (hamming_distance == 0.0) begin
                $display("PERFECT: Exact key recovered!");
            end else if (hamming_distance <= 8.0) begin
                $display("EXCELLENT: Very close to true key");
            end else if (hamming_distance <= 16.0) begin
                $display("GOOD: Reasonably close");
            end else if (correct_bytes > 0) begin
                $display("PARTIAL: Some bytes recovered");
            end else begin
                $display("POOR: No bytes recovered");
            end
            
            if (best_correlation > 0.9) begin
                $display("Attack achieved high correlation - likely successful");
            end else if (best_correlation > 0.5) begin
                $display("Attack achieved moderate correlation");
            end else begin
                $display("Attack achieved low correlation");
            end
        end
    endtask
    
    function real calculate_hamming_distance;
        input [63:0] a, b;
        integer diff, count;
        begin
            diff = a ^ b;
            count = 0;
            while (diff != 0) begin
                if (diff & 1) count = count + 1;
                diff = diff >> 1;
            end
            calculate_hamming_distance = $itor(count);
        end
    endfunction
    
    function integer count_correct_bytes;
        input [63:0] a, b;
        integer correct;
        begin
            correct = 0;
            if (a[63:56] == b[63:56]) correct = correct + 1;
            if (a[55:48] == b[55:48]) correct = correct + 1;
            if (a[47:40] == b[47:40]) correct = correct + 1;
            if (a[39:32] == b[39:32]) correct = correct + 1;
            if (a[31:24] == b[31:24]) correct = correct + 1;
            if (a[23:16] == b[23:16]) correct = correct + 1;
            if (a[15:8] == b[15:8]) correct = correct + 1;
            if (a[7:0] == b[7:0]) correct = correct + 1;
            count_correct_bytes = correct;
        end
    endfunction
    
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
                1: key[15:8] = value[7:0];
                0: key[7:0] = value[7:0];
            endcase
        end
    endtask
    
    task complete_system_reset;
        begin
            cs = 0;
            we = 0;
            addr = 0;
            wdata = 0;
            trng_a_in = 0;
            trng_d_in = 0;
            dcr = 0;
            
            repeat (20) @(posedge clk);
            
            rst_n = 0;
            repeat (5) @(posedge clk);
            rst_n = 1;
            repeat (10) @(posedge clk);
        end
    endtask

endmodule