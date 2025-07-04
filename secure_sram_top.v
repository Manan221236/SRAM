// `timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.06.2025 14:46:19
// Design Name: 
// Module Name: secure_sram_top
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////

module secure_sram_top #(
    parameter ADDR_WIDTH = 14,
    parameter DATA_WIDTH = 52,
    parameter TRNG_A_WIDTH = 64,  // TRNG for address mapping
    parameter TRNG_D_WIDTH = 32,  // TRNG for data mapping
    parameter MEMORY_DEPTH = 2**ADDR_WIDTH
)(
    input wire clk,
    input wire rst_n,
    input wire dcr,                    // Data Corruption Request
    input wire cs,                     // Chip Select
    input wire we,                     // Write Enable
    input wire [ADDR_WIDTH-1:0] addr,
    input wire [DATA_WIDTH-1:0] wdata,
    input wire [TRNG_A_WIDTH-1:0] trng_a_in,
    input wire [TRNG_D_WIDTH-1:0] trng_d_in,
    output reg [DATA_WIDTH-1:0] rdata,
    output wire ready
);

    // Internal signals
    wire [ADDR_WIDTH-1:0] mapped_addr;
    wire [DATA_WIDTH-1:0] mapped_wdata, demapped_rdata;
    wire [DATA_WIDTH-1:0] sram_rdata;
    wire [TRNG_A_WIDTH-1:0] trng_a_reg;
    wire [TRNG_D_WIDTH-1:0] trng_d_reg;
    
    reg [1:0] pipeline_stage;
    reg cs_reg, we_reg;
    
    // TRNG Register Management
    trng_register_manager #(
        .TRNG_A_WIDTH(TRNG_A_WIDTH),
        .TRNG_D_WIDTH(TRNG_D_WIDTH)
    ) trng_mgr (
        .clk(clk),
        .rst_n(rst_n),
        .dcr(dcr),
        .trng_a_in(trng_a_in),
        .trng_d_in(trng_d_in),
        .trng_a_out(trng_a_reg),
        .trng_d_out(trng_d_reg)
    );
    
    // Address Mapping
    address_mapper #(
        .ADDR_WIDTH(ADDR_WIDTH),
        .TRNG_WIDTH(TRNG_A_WIDTH)
    ) addr_map (
        .addr_in(addr),
        .trng(trng_a_reg),
        .addr_out(mapped_addr)
    );
    
    // Data Mapping (for write)
    data_mapper #(
        .DATA_WIDTH(DATA_WIDTH),
        .TRNG_WIDTH(TRNG_D_WIDTH)
    ) data_map (
        .data_in(wdata),
        .trng(trng_d_reg),
        .data_out(mapped_wdata)
    );
    
    // Data De-mapping (for read)
    data_mapper #(
        .DATA_WIDTH(DATA_WIDTH),
        .TRNG_WIDTH(TRNG_D_WIDTH)
    ) data_demap (
        .data_in(sram_rdata),
        .trng(trng_d_reg),
        .data_out(demapped_rdata)
    );
    
    // SRAM Array
    sram_array #(
        .ADDR_WIDTH(ADDR_WIDTH),
        .DATA_WIDTH(DATA_WIDTH),
        .MEMORY_DEPTH(MEMORY_DEPTH)
    ) sram_core (
        .clk(clk),
        .cs(cs_reg),
        .we(we_reg),
        .addr(mapped_addr),
        .wdata(mapped_wdata),
        .rdata(sram_rdata)
    );
    
    // Pipeline control for timing
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            pipeline_stage <= 2'b00;
            cs_reg <= 1'b0;
            we_reg <= 1'b0;
            rdata <= {DATA_WIDTH{1'b0}};
        end else begin
            // Stage 1: Register inputs and perform mapping
            cs_reg <= cs;
            we_reg <= we;
            
            // Stage 2: SRAM access and de-mapping
            if (cs_reg && !we_reg) begin
                rdata <= demapped_rdata;
            end
        end
    end
    
    assign ready = 1'b1; // Simple ready signal

endmodule

