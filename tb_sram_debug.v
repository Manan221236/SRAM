`timescale 1ns / 1ps

module tb_sram_debug;
    
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
    
    reg [7:0] test_pixels [0:15];
    reg [7:0] read_pixels [0:15];
    
    reg [63:0] test_trng_a;
    reg [31:0] test_trng_d;
    
    integer i;
    integer matches;
    
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
        $display("SRAM SYSTEM DIAGNOSTIC");
        $display("======================");
        
        // Initialize
        rst_n = 0;
        dcr = 0;
        cs = 0;
        we = 0;
        addr = 0;
        wdata = 0;
        trng_a_in = 0;
        trng_d_in = 0;
        
        test_trng_a = 64'hDEADBEEFCAFEBABE;
        test_trng_d = 32'h12345678;
        
        repeat (10) @(posedge clk);
        rst_n = 1;
        repeat (5) @(posedge clk);
        
        // Test 1: Basic functionality
        test_basic_operation;
        
        // Test 2: Key dependency
        test_key_dependency;
        
        // Test 3: DCR impact
        test_dcr_behavior;
        
        $finish;
    end

    task test_basic_operation;
        begin
            $display("\n=== TEST 1: Basic Operation ===");
            
            // Create simple test pattern
            for (i = 0; i < 16; i = i + 1) begin
                test_pixels[i] = i * 16 + 5;  // 5, 21, 37, 53, ...
            end
            
            $display("Original pixels: %0d %0d %0d %0d %0d %0d %0d %0d", 
                     test_pixels[0], test_pixels[1], test_pixels[2], test_pixels[3],
                     test_pixels[4], test_pixels[5], test_pixels[6], test_pixels[7]);
            
            // Reset system
            complete_system_reset;
            
            // Configure keys
            trng_a_in = test_trng_a;
            trng_d_in = test_trng_d;
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (5) @(posedge clk);
            
            // Write test data
            $display("Writing data...");
            for (i = 0; i < 16; i = i + 1) begin
                @(posedge clk);
                cs = 1;
                we = 1;
                addr = i;
                wdata = {43'd0, test_pixels[i]};
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                cs = 0;
                we = 0;
            end
            
            $display("Data written. Now reading back...");
            
            // Read back data with SAME keys
            for (i = 0; i < 16; i = i + 1) begin
                @(posedge clk);
                cs = 1;
                we = 0;
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                read_pixels[i] = rdata[7:0];
                cs = 0;
            end
            
            $display("Read pixels:     %0d %0d %0d %0d %0d %0d %0d %0d", 
                     read_pixels[0], read_pixels[1], read_pixels[2], read_pixels[3],
                     read_pixels[4], read_pixels[5], read_pixels[6], read_pixels[7]);
            
            // Check match
            matches = 0;
            for (i = 0; i < 16; i = i + 1) begin
                if (test_pixels[i] == read_pixels[i]) matches = matches + 1;
            end
            
            $display("Result: %0d/16 pixels match", matches);
            
            if (matches == 16) begin
                $display("PASS: Perfect round-trip with same keys");
            end else if (matches > 0) begin
                $display("PARTIAL: Some pixels match - encryption may be weak");
            end else begin
                $display("FAIL: No pixels match - system broken or very strong encryption");
            end
        end
    endtask

    task test_key_dependency;
        begin
            $display("\n=== TEST 2: Key Dependency ===");
            
            // Use same test data, but read with DIFFERENT keys
            complete_system_reset;
            
            // Try with different TRNG_A
            trng_a_in = 64'h0000000000000000;  // Different from test_trng_a
            trng_d_in = test_trng_d;           // Same TRNG_D
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (5) @(posedge clk);
            
            $display("Reading with wrong TRNG_A (zeros)...");
            for (i = 0; i < 16; i = i + 1) begin
                @(posedge clk);
                cs = 1;
                we = 0;
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                read_pixels[i] = rdata[7:0];
                cs = 0;
            end
            
            $display("Wrong key pixels: %0d %0d %0d %0d %0d %0d %0d %0d", 
                     read_pixels[0], read_pixels[1], read_pixels[2], read_pixels[3],
                     read_pixels[4], read_pixels[5], read_pixels[6], read_pixels[7]);
            
            matches = 0;
            for (i = 0; i < 16; i = i + 1) begin
                if (test_pixels[i] == read_pixels[i]) matches = matches + 1;
            end
            
            $display("Wrong TRNG_A result: %0d/16 pixels match", matches);
            
            // Try with different TRNG_D
            complete_system_reset;
            trng_a_in = test_trng_a;           // Correct TRNG_A
            trng_d_in = 32'h87654321;          // Different TRNG_D
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (5) @(posedge clk);
            
            $display("Reading with wrong TRNG_D...");
            for (i = 0; i < 16; i = i + 1) begin
                @(posedge clk);
                cs = 1;
                we = 0;
                addr = i;
                @(posedge clk);
                while (!ready) @(posedge clk);
                @(posedge clk);
                read_pixels[i] = rdata[7:0];
                cs = 0;
            end
            
            $display("Wrong TRNG_D pixels: %0d %0d %0d %0d %0d %0d %0d %0d", 
                     read_pixels[0], read_pixels[1], read_pixels[2], read_pixels[3],
                     read_pixels[4], read_pixels[5], read_pixels[6], read_pixels[7]);
            
            matches = 0;
            for (i = 0; i < 16; i = i + 1) begin
                if (test_pixels[i] == read_pixels[i]) matches = matches + 1;
            end
            
            $display("Wrong TRNG_D result: %0d/16 pixels match", matches);
            
            if (matches == 0) begin
                $display("GOOD: Keys matter - encryption is working");
            end else begin
                $display("CONCERN: Wrong keys still give some matches");
            end
        end
    endtask

    task test_dcr_behavior;
        begin
            $display("\n=== TEST 3: DCR Impact ===");
            
            // Test what happens without DCR pulse
            complete_system_reset;
            
            trng_a_in = test_trng_a;
            trng_d_in = test_trng_d;
            // Skip DCR pulse
            repeat (5) @(posedge clk);
            
            $display("Writing without DCR pulse...");
            @(posedge clk);
            cs = 1;
            we = 1;
            addr = 0;
            wdata = {43'd0, 8'hAA};  // Test value
            @(posedge clk);
            while (!ready) @(posedge clk);
            @(posedge clk);
            cs = 0;
            we = 0;
            
            $display("Reading back...");
            @(posedge clk);
            cs = 1;
            we = 0;
            addr = 0;
            @(posedge clk);
            while (!ready) @(posedge clk);
            @(posedge clk);
            $display("Without DCR: wrote 0xAA, read 0x%02X", rdata[7:0]);
            cs = 0;
            
            // Now test WITH DCR pulse
            complete_system_reset;
            
            trng_a_in = test_trng_a;
            trng_d_in = test_trng_d;
            dcr = 1;
            repeat (10) @(posedge clk);
            dcr = 0;
            repeat (5) @(posedge clk);
            
            @(posedge clk);
            cs = 1;
            we = 1;
            addr = 0;
            wdata = {43'd0, 8'hAA};  // Same test value
            @(posedge clk);
            while (!ready) @(posedge clk);
            @(posedge clk);
            cs = 0;
            we = 0;
            
            @(posedge clk);
            cs = 1;
            we = 0;
            addr = 0;
            @(posedge clk);
            while (!ready) @(posedge clk);
            @(posedge clk);
            $display("With DCR:    wrote 0xAA, read 0x%02X", rdata[7:0]);
            cs = 0;
            
            if (rdata[7:0] == 8'hAA) begin
                $display("DCR enables transparent mode");
            end else begin
                $display("DCR enables encryption mode");
            end
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