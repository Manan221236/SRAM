`timescale 1ns / 1ps

module tb_attack_multi_start_search;
    
    parameter NUM_ITERATIONS = 10;
    parameter NUM_RANDOM_STARTS = 8;
    parameter STOP_WHEN_PERFECT = 1;
    parameter DUMP_EVERY_BYTE = 0;
    parameter DUMP_EVERY_ITER = 1;
    parameter ENABLE_DIAGNOSTICS = 1;
    parameter ENABLE_KEY_VERIFICATION = 1;
    parameter RANDOM_SEED = 42;
    
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
    
    reg [7:0] plain_pix [0:16383];
    reg [7:0] decoded_pix [0:16383];
    reg [7:0] temp_decoded [0:16383];
    
    reg [63:0] true_trng_a;
    reg [31:0] true_trng_d;
    reg [31:0] bad_trng_d;
    reg [63:0] found_trng_a;
    reg [63:0] global_best_key;
    real global_best_score;
    
    reg [63:0] start_keys [0:NUM_RANDOM_STARTS-1];
    reg [63:0] start_results [0:NUM_RANDOM_STARTS-1];
    real start_scores [0:NUM_RANDOM_STARTS-1];
    integer best_start_idx;
    
    reg [63:0] iteration_keys [0:NUM_ITERATIONS-1];
    real iteration_composite [0:NUM_ITERATIONS-1];
    
    integer i, j, iter_idx, byte_idx, val, start_idx;
    integer fh;
    reg [63:0] cand;
    real cur_pearson, cur_spearman, cur_cosine, cur_mse, cur_composite;
    real best_composite, prev_score;
    integer best_val;
    reg [8*50:1] fname;
    
    real test_p, test_s, test_c, test_m, test_comp;
    real pre_search_score, post_search_score;
    integer consistency_errors;
    
    real true_key_score, wrong_key_score, zero_key_score;
    
    reg [31:0] rng_state;
    
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
        $display("Starting MULTI-START RANDOM SEARCH Attack Testbench");
        
        rst_n = 0;
        dcr = 0;
        cs = 0;
        we = 0;
        addr = 0;
        wdata = 0;
        trng_a_in = 0;
        trng_d_in = 0;
        consistency_errors = 0;
        global_best_score = -1000.0;
        global_best_key = 64'h0;
        best_start_idx = 0;
        rng_state = RANDOM_SEED;
        
        load_realistic_test_image;
        
        true_trng_a = 64'hDEADBEEFCAFEBABE;
        true_trng_d = 32'h12345678;
        bad_trng_d = 32'h87654321;
        
        repeat (10) @(posedge clk);
        rst_n = 1;
        repeat (5) @(posedge clk);
        
        encode_reference_image;
        decode_reference_image;
        generate_comparison_case;
        
        if (ENABLE_KEY_VERIFICATION) begin
            perform_key_verification;
        end
        
        $display("Starting multi-start search with %0d starting points", NUM_RANDOM_STARTS);
        
        generate_random_starting_points;
        
        for (start_idx = 0; start_idx < NUM_RANDOM_STARTS; start_idx = start_idx + 1) begin
            $display("Starting point %0d: %h", start_idx, start_keys[start_idx]);
            
            found_trng_a = start_keys[start_idx];
            
            for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
                $display("Start %0d - Iteration %0d", start_idx, iter_idx);
                
                for (byte_idx = 7; byte_idx >= 0; byte_idx = byte_idx - 1) begin
                    $display("Searching byte %0d", byte_idx);
                    
                    if (ENABLE_DIAGNOSTICS) begin
                        corrected_score_candidate(found_trng_a, test_p, test_s, test_c, test_m, pre_search_score);
                    end
                    
                    best_composite = -1000.0;
                    best_val = get_current_byte(found_trng_a, byte_idx);
                    prev_score = pre_search_score;
                    
                    for (val = 0; val < 256; val = val + 1) begin
                        cand = found_trng_a;
                        set_byte_value(cand, byte_idx, val);
                        
                        corrected_score_candidate(cand, cur_pearson, cur_spearman, 
                                               cur_cosine, cur_mse, cur_composite);
                        
                        if (ENABLE_DIAGNOSTICS && val == get_current_byte(found_trng_a, byte_idx)) begin
                            if (cur_composite > pre_search_score) begin
                                if ((cur_composite - pre_search_score) > 0.000001) begin
                                    consistency_errors = consistency_errors + 1;
                                end
                            end else begin
                                if ((pre_search_score - cur_composite) > 0.000001) begin
                                    consistency_errors = consistency_errors + 1;
                                end
                            end
                        end
                        
                        if (cur_composite > best_composite) begin
                            best_composite = cur_composite;
                            best_val = val;
                        end
                        
                        if (STOP_WHEN_PERFECT && cur_composite >= 0.95) begin
                            val = 256;
                        end
                    end
                    
                    if (best_composite > prev_score) begin
                        if ((best_composite - prev_score) > 0.000001) begin
                            set_byte_value(found_trng_a, byte_idx, best_val);
                            $display("Improved byte %0d to 0x%02X (score: %0.6f)", 
                                     byte_idx, best_val[7:0], best_composite);
                        end
                    end
                    
                    if (DUMP_EVERY_BYTE) dump_byte_snapshot(start_idx, iter_idx, byte_idx);
                end
                
                if (start_idx == 0) begin
                    iteration_keys[iter_idx] = found_trng_a;
                    corrected_score_candidate(found_trng_a, test_p, test_s, test_c, test_m, iteration_composite[iter_idx]);
                end
                
                if (DUMP_EVERY_ITER) dump_iteration_snapshot(start_idx, iter_idx);
            end
            
            start_results[start_idx] = found_trng_a;
            corrected_score_candidate(found_trng_a, test_p, test_s, test_c, test_m, start_scores[start_idx]);
            
            $display("Start %0d final: %h (score: %0.6f)", 
                     start_idx, start_results[start_idx], start_scores[start_idx]);
            
            if (start_scores[start_idx] > global_best_score) begin
                global_best_score = start_scores[start_idx];
                global_best_key = start_results[start_idx];
                best_start_idx = start_idx;
                $display("New global best from start %0d: %0.6f", start_idx, global_best_score);
            end
        end
        
        $display("Multi-start search complete");
        $display("Best result from start %0d: %0.6f", best_start_idx, global_best_score);
        $display("Best key: %h", global_best_key);
        
        if (true_key_score <= global_best_score) begin
            $display("ERROR: Attack found key better than ground truth!");
        end else begin
            $display("PASS: True key has highest score");
        end
        
        found_trng_a = global_best_key;
        for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
            iteration_keys[iter_idx] = global_best_key;
            iteration_composite[iter_idx] = global_best_score;
        end
        
        dump_final_comparison;
        generate_summary_file;
        $finish;
    end
    
    function [31:0] next_random;
        input [31:0] current_state;
        begin
            next_random = {current_state[30:0], current_state[31] ^ current_state[21] ^ current_state[1] ^ current_state[0]};
        end
    endfunction
    
    task generate_random_starting_points;
        integer temp_rng;
        reg [31:0] low_part, high_part;
        begin
            temp_rng = rng_state;
            
            for (start_idx = 0; start_idx < NUM_RANDOM_STARTS; start_idx = start_idx + 1) begin
                temp_rng = next_random(temp_rng);
                low_part = temp_rng;
                temp_rng = next_random(temp_rng);
                high_part = temp_rng;
                start_keys[start_idx] = {high_part, low_part};
            end
            
            if (NUM_RANDOM_STARTS >= 4) begin
                start_keys[NUM_RANDOM_STARTS-4] = 64'h0000000000000000;
                start_keys[NUM_RANDOM_STARTS-3] = 64'hFFFFFFFFFFFFFFFF;
                start_keys[NUM_RANDOM_STARTS-2] = 64'hAAAAAAAAAAAAAAAA;
                start_keys[NUM_RANDOM_STARTS-1] = 64'h5555555555555555;
            end
        end
    endtask
    
    task perform_key_verification;
        reg [63:0] test_key;
        begin
            $display("KEY VERIFICATION - Testing Known Keys");
            
            test_key = true_trng_a;
            corrected_score_candidate(test_key, test_p, test_s, test_c, test_m, true_key_score);
            $display("TRUE key: %h scores %0.6f", test_key, true_key_score);
            
            test_key = 64'h0;
            corrected_score_candidate(test_key, test_p, test_s, test_c, test_m, zero_key_score);
            $display("ZERO key: scores %0.6f", zero_key_score);
            
            test_key = 64'hFFFFFFFFFFFFFFFF;
            corrected_score_candidate(test_key, test_p, test_s, test_c, test_m, wrong_key_score);
            $display("WRONG key: scores %0.6f", wrong_key_score);
            
            if (true_key_score > zero_key_score && true_key_score > wrong_key_score) begin
                $display("PASS: True key has highest score");
            end else begin
                $display("FAIL: True key does not have highest score");
            end
        end
    endtask
    
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
    
    task corrected_score_candidate;
        input [63:0] key;
        output real pearson_score;
        output real spearman_score;
        output real cosine_score;
        output real mse_score;
        output real composite_score;
        
        integer sum_x, sum_y, sum_xy, sum_x2, sum_y2;
        integer n, numerator, denom_x, denom_y;
        real denom_real, dot_product, norm_x, norm_y, mse_raw;
        integer temp_val;
        real raw_pearson, raw_cosine, raw_mse_sim;
        
        begin
            complete_system_reset;
            
            repeat (20) @(posedge clk);
            
            trng_a_in = key;
            trng_d_in = bad_trng_d;
            
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (10) @(posedge clk);
            
            for (i = 0; i < 16384; i = i + 1) begin
                repeat (2) @(posedge clk);
                cs = 1;
                we = 0;
                addr = i;
                repeat (2) @(posedge clk);
                while (!ready) @(posedge clk);
                repeat (2) @(posedge clk);
                temp_decoded[i] = rdata[7:0];
                cs = 0;
                repeat (2) @(posedge clk);
            end
            
            cs = 0;
            we = 0;
            addr = 0;
            repeat (10) @(posedge clk);
            
            for (i = 0; i < 16384; i = i + 1) begin
                decoded_pix[i] = temp_decoded[i];
            end
            
            n = 16384;
            
            sum_x = 0; sum_y = 0; sum_xy = 0; sum_x2 = 0; sum_y2 = 0;
            
            for (i = 0; i < n; i = i + 1) begin
                sum_x = sum_x + plain_pix[i];
                sum_y = sum_y + decoded_pix[i];
                sum_xy = sum_xy + (plain_pix[i] * decoded_pix[i]);
                sum_x2 = sum_x2 + (plain_pix[i] * plain_pix[i]);
                sum_y2 = sum_y2 + (decoded_pix[i] * decoded_pix[i]);
            end
            
            numerator = (n * sum_xy) - (sum_x * sum_y);
            denom_x = (n * sum_x2) - (sum_x * sum_x);
            denom_y = (n * sum_y2) - (sum_y * sum_y);
            
            if (denom_x <= 0 || denom_y <= 0) begin
                raw_pearson = 0.0;
            end else begin
                denom_real = $sqrt($itor(denom_x) * $itor(denom_y));
                raw_pearson = $itor(numerator) / denom_real;
            end
            
            if (raw_pearson > 1.0) raw_pearson = 1.0;
            if (raw_pearson < -1.0) raw_pearson = -1.0;
            
            spearman_score = raw_pearson;
            
            dot_product = 0.0;
            norm_x = 0.0;
            norm_y = 0.0;
            
            for (i = 0; i < n; i = i + 1) begin
                dot_product = dot_product + ($itor(plain_pix[i]) * $itor(decoded_pix[i]));
                norm_x = norm_x + ($itor(plain_pix[i]) * $itor(plain_pix[i]));
                norm_y = norm_y + ($itor(decoded_pix[i]) * $itor(decoded_pix[i]));
            end
            
            norm_x = $sqrt(norm_x);
            norm_y = $sqrt(norm_y);
            
            if (norm_x == 0.0 || norm_y == 0.0) begin
                raw_cosine = 0.0;
            end else begin
                raw_cosine = dot_product / (norm_x * norm_y);
            end
            
            if (raw_cosine > 1.0) raw_cosine = 1.0;
            if (raw_cosine < -1.0) raw_cosine = -1.0;
            
            mse_raw = 0.0;
            for (i = 0; i < n; i = i + 1) begin
                temp_val = plain_pix[i] - decoded_pix[i];
                mse_raw = mse_raw + ($itor(temp_val) * $itor(temp_val));
            end
            mse_raw = mse_raw / $itor(n);
            
            raw_mse_sim = 1.0 - (mse_raw / 65025.0);
            if (raw_mse_sim < 0.0) raw_mse_sim = 0.0;
            if (raw_mse_sim > 1.0) raw_mse_sim = 1.0;
            
            pearson_score = (raw_pearson + 1.0) / 2.0;
            spearman_score = (raw_pearson + 1.0) / 2.0;
            cosine_score = (raw_cosine + 1.0) / 2.0;
            mse_score = raw_mse_sim;
            
            if (pearson_score < 0.0) pearson_score = 0.0;
            if (pearson_score > 1.0) pearson_score = 1.0;
            if (cosine_score < 0.0) cosine_score = 0.0;
            if (cosine_score > 1.0) cosine_score = 1.0;
            if (mse_score < 0.0) mse_score = 0.0;
            if (mse_score > 1.0) mse_score = 1.0;
            
            composite_score = 0.50 * pearson_score + 0.20 * cosine_score + 0.30 * mse_score;
            
            if (composite_score < 0.0) composite_score = 0.0;
            if (composite_score > 1.0) composite_score = 1.0;
            
            complete_system_reset;
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
    
    task load_realistic_test_image;
        integer seed;
        integer temp_val;
        begin
            seed = 42;
            
            for (i = 0; i < 16384; i = i + 1) begin
                temp_val = ((i % 256) ^ ((i/128) % 256) ^ ((i*3) % 256));
                plain_pix[i] = temp_val[7:0];
                
                if (i % 1000 == 0) begin
                    plain_pix[i] = plain_pix[i] ^ 8'h55;
                end
            end
            $display("Loaded test image");
        end
    endtask
    
    task encode_reference_image;
        begin
            complete_system_reset;
            
            trng_a_in = true_trng_a;
            trng_d_in = true_trng_d;
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (10) @(posedge clk);

            for (i = 0; i < 16384; i = i + 1) begin
                @(posedge clk);
                cs = 1;
                we = 1;
                addr = i;
                wdata = {43'd0, plain_pix[i]};
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                cs = 0;
            end
            $display("Encoded reference image");
            complete_system_reset;
        end
    endtask
    
    task decode_reference_image;
        begin
            complete_system_reset;
            
            trng_a_in = true_trng_a;
            trng_d_in = true_trng_d;
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (10) @(posedge clk);

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
            $display("Generated reference decode file");
            complete_system_reset;
        end
    endtask
    
    task generate_comparison_case;
        begin
            complete_system_reset;
            
            trng_a_in = true_trng_a;
            trng_d_in = bad_trng_d;
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (10) @(posedge clk);

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
            $display("Generated comparison case file");
            complete_system_reset;
        end
    endtask
    
    task dump_byte_snapshot;
        input integer start_num;
        input integer iter_num;
        input integer byte_num;
        begin
            $sformat(fname, "byte_s%0d_i%0d_b%0d.txt", start_num, iter_num, byte_num);
            fh = $fopen(fname,"w");
            for (i = 0; i < 16384; i = i + 1) begin
                $fwrite(fh,"%013X\n", {43'd0, decoded_pix[i]});
            end
            $fclose(fh);
        end
    endtask
    
    task dump_iteration_snapshot;
        input integer start_num;
        input integer iter_num;
        begin
            $sformat(fname, "iter_s%0d_i%0d.txt", start_num, iter_num);
            fh = $fopen(fname,"w");
            for (i = 0; i < 16384; i = i + 1) begin
                $fwrite(fh,"%013X\n", {43'd0, decoded_pix[i]});
            end
            $fclose(fh);
        end
    endtask
    
    task dump_final_comparison;
        begin
            for (iter_idx = 0; iter_idx < NUM_ITERATIONS; iter_idx = iter_idx + 1) begin
                complete_system_reset;
                
                trng_a_in = global_best_key;
                trng_d_in = bad_trng_d;
                dcr = 1;
                repeat (10) @(posedge clk);
                dcr = 0;
                repeat (10) @(posedge clk);

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
                    $fwrite(fh,"%013X\n", rdata);
                    cs = 0;
                end
                $fclose(fh);
            end
            $display("Generated final comparison files");
            complete_system_reset;
        end
    endtask
    
    task generate_summary_file;
        begin
            fh = $fopen("multi_metric_summary.txt","w");
            $fwrite(fh,"Multi-Start Random Search Summary\n");
            $fwrite(fh,"=================================\n\n");
            $fwrite(fh,"Ground Truth: %h\n", true_trng_a);
            $fwrite(fh,"Random starts: %0d\n", NUM_RANDOM_STARTS);
            $fwrite(fh,"Iterations: %0d\n", NUM_ITERATIONS);
            $fwrite(fh,"Errors: %0d\n\n", consistency_errors);
            
            $fwrite(fh,"Key Verification:\n");
            $fwrite(fh,"True key: %0.6f\n", true_key_score);
            $fwrite(fh,"Zero key: %0.6f\n", zero_key_score);
            $fwrite(fh,"Wrong key: %0.6f\n\n", wrong_key_score);
            
            $fwrite(fh,"Results:\n");
            for (start_idx = 0; start_idx < NUM_RANDOM_STARTS; start_idx = start_idx + 1) begin
                if (start_idx == best_start_idx) begin
                    $fwrite(fh,"Start %0d: %h (%0.6f) BEST\n", 
                            start_idx, start_results[start_idx], start_scores[start_idx]);
                end else begin
                    $fwrite(fh,"Start %0d: %h (%0.6f)\n", 
                            start_idx, start_results[start_idx], start_scores[start_idx]);
                end
            end
            
            $fwrite(fh,"\nBest: %h (%0.6f)\n", global_best_key, global_best_score);
            
            if (true_key_score > global_best_score) begin
                $fwrite(fh,"Status: PASS - True key highest\n");
            end else begin
                $fwrite(fh,"Status: FAIL - Attack too good\n");
            end
            
            $fclose(fh);
            $display("Generated summary file");
        end
    endtask

endmodule