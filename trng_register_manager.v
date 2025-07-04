`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.06.2025 14:46:19
// Design Name: 
// Module Name: trng_register_manager
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


module trng_register_manager #(
    parameter TRNG_A_WIDTH = 64,
    parameter TRNG_D_WIDTH = 32
)(
    input wire clk,
    input wire rst_n,
    input wire dcr,
    input wire [TRNG_A_WIDTH-1:0] trng_a_in,
    input wire [TRNG_D_WIDTH-1:0] trng_d_in,
    output reg [TRNG_A_WIDTH-1:0] trng_a_out,
    output reg [TRNG_D_WIDTH-1:0] trng_d_out
);

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            trng_a_out <= {TRNG_A_WIDTH{1'b0}};
            trng_d_out <= {TRNG_D_WIDTH{1'b0}};
        end else if (dcr) begin
            // Load new TRNGs on DCR event (single cycle)
            trng_a_out <= trng_a_in;
            trng_d_out <= trng_d_in;
        end
    end

endmodule
