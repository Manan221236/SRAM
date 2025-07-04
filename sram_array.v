`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.06.2025 14:46:19
// Design Name: 
// Module Name: sram_array
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


module sram_array #(
    parameter ADDR_WIDTH = 14,
    parameter DATA_WIDTH = 52,
    parameter MEMORY_DEPTH = 2**ADDR_WIDTH
)(
    input wire clk,
    input wire cs,
    input wire we,
    input wire [ADDR_WIDTH-1:0] addr,
    input wire [DATA_WIDTH-1:0] wdata,
    output reg [DATA_WIDTH-1:0] rdata
);

    // Memory array
    reg [DATA_WIDTH-1:0] memory [0:MEMORY_DEPTH-1];
    initial begin
        for (integer i = 0; i < MEMORY_DEPTH; i = i + 1)
            memory[i] = {DATA_WIDTH{1'b0}};
    end
    always @(posedge clk) begin
        if (cs) begin
            if (we) begin
                // Write operation
                memory[addr] <= wdata;
            end else begin
                // Read operation
                rdata <= memory[addr];
            end
        end
    end

endmodule
