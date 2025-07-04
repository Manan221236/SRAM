// `timescale 1ns/1ps
/*  =======================  USER SWITCHES  ============================
 *  Edit if you want to change behaviour (no SV needed; plain defines)
 *  ------------------------------------------------------------------*/
`define STOP_WHEN_PERFECT   1   // 1 = break inner loop when correlation >= 0.95
`define DUMP_EVERY_BYTE     1   // 0 = only dump final snapshot
`define CORRELATION_THRESH  950 // 0.95 * 1000 (scaled integer)
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
    parameter CORRELATION_THRESH = `CORRELATION_THRESH;

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
    
    /* ------------ locals ------------------------------------------- */
    integer i, val, best_val, best_corr, cur_corr;
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

        /* ---------- iterative brute-force per byte (MSB→LSB) -------- */
        found_trng_a = 64'h0;

        for (byte_idx = 7; byte_idx >= 0; byte_idx = byte_idx - 1) begin
            $display("\n=== Searching byte %0d ===", byte_idx);
            best_corr = -1000; 
            best_val = 8'h00;  // correlation * 1000 (can be negative)

            for (val = 0; val < 256; val = val + 1) begin : SEARCH
                cand = found_trng_a;
                // Traditional Verilog bit selection
                case (byte_idx)
                    7: cand[63:56] = val[7:0];
                    6: cand[55:48] = val[7:0];
                    5: cand[47:40] = val[7:0];
                    4: cand[39:32] = val[7:0];
                    3: cand[31:24] = val[7:0];
                    2: cand[23:16] = val[7:0];
                    1: cand[15:8]  = val[7:0];
                    0: cand[7:0]   = val[7:0];
                endcase
                
                pearson_score_candidate(cand, cur_corr);

                if (cur_corr > best_corr) begin
                    best_corr = cur_corr;
                    best_val  = val;
                end
                if (STOP_WHEN_PERFECT && cur_corr >= CORRELATION_THRESH)
                    disable SEARCH;               // early exit
            end

            // Traditional Verilog bit assignment
            case (byte_idx)
                7: found_trng_a[63:56] = best_val[7:0];
                6: found_trng_a[55:48] = best_val[7:0];
                5: found_trng_a[47:40] = best_val[7:0];
                4: found_trng_a[39:32] = best_val[7:0];
                3: found_trng_a[31:24] = best_val[7:0];
                2: found_trng_a[23:16] = best_val[7:0];
                1: found_trng_a[15:8]  = best_val[7:0];
                0: found_trng_a[7:0]   = best_val[7:0];
            endcase
            
            $display("byte %0d chosen = 0x%02X  (correlation %0d.%03d)",
                     byte_idx, best_val[7:0], best_corr/1000, 
                     (best_corr >= 0) ? (best_corr % 1000) : ((-best_corr) % 1000));

            if (DUMP_EVERY_BYTE) dump_snapshot(byte_idx);
        end

        $display("\nRecovered TRNG_A = %h", found_trng_a);
        $display("Ground-truth      = %h", true_trng_a);

        dump_snapshot(8);   // final image with full recovered key
        $finish;
    end

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
        output integer correlation_x1000;
        
        // Correlation calculation variables
        integer sum_x, sum_y, sum_xy, sum_x2, sum_y2;
        integer n;
        integer numerator, denom_x, denom_y;
        real correlation_real;
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
            
            // Calculate correlation using the formula:
            // r = (n*Σxy - Σx*Σy) / sqrt((n*Σx² - (Σx)²) * (n*Σy² - (Σy)²))
            
            numerator = (n * sum_xy) - (sum_x * sum_y);
            denom_x   = (n * sum_x2) - (sum_x * sum_x);
            denom_y   = (n * sum_y2) - (sum_y * sum_y);
            
            // Handle edge cases
            if (denom_x <= 0 || denom_y <= 0) begin
                correlation_x1000 = -1000;  // Invalid correlation
            end else begin
                // Calculate denominator using real arithmetic
                denom_real = $sqrt($itor(denom_x) * $itor(denom_y));
                correlation_real = $itor(numerator) / denom_real;
                correlation_x1000 = $rtoi(correlation_real * 1000.0);
                
                // Clamp to valid range [-1000, 1000]
                if (correlation_x1000 > 1000) correlation_x1000 = 1000;
                if (correlation_x1000 < -1000) correlation_x1000 = -1000;
            end
        end
    endtask

    /* =================================================================
     *  TASK: dump_snapshot
     *        save whole SRAM view with current found_trng_a
     *        idx 7..0 = per-byte, idx 8 = final
     * =================================================================*/
    task dump_snapshot;
        input integer idx;
        begin
            trng_a_in = found_trng_a;
            trng_d_in = bad_trng_d;
            dcr = 1; 
            repeat (3) @(posedge clk); 
            dcr = 0; 
            @(posedge clk);

            if (idx == 8)
                $sformat(fname, "recover_final.txt");
            else
                $sformat(fname, "recover_byte%0d.txt", idx);
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

endmodule