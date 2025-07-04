`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.06.2025 14:46:19
// Design Name: 
// Module Name: simple_trng
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


module simple_trng #(
    parameter WIDTH = 32,
    parameter SEED = 32'hACE1_BEEF
)(
    input wire clk,
    input wire rst_n,
    input wire enable,
    output reg [WIDTH-1:0] trng_out
);

    reg [WIDTH-1:0] lfsr;
    wire feedback;
    
    // Simple LFSR feedback for 32-bit (polynomial: x^32 + x^22 + x^2 + x^1 + 1)
    assign feedback = lfsr[31] ^ lfsr[21] ^ lfsr[1] ^ lfsr[0];
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            lfsr <= SEED;
            trng_out <= SEED;
        end else if (enable) begin
            lfsr <= {lfsr[WIDTH-2:0], feedback};
            trng_out <= lfsr;
        end
    end

endmodule
