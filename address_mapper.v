`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.06.2025 14:46:19
// Design Name: 
// Module Name: address_mapper
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
module address_mapper #(
    parameter ADDR_WIDTH = 14,
    parameter TRNG_WIDTH = 64
)(
    input wire [ADDR_WIDTH-1:0] addr_in,
    input wire [TRNG_WIDTH-1:0] trng,
    output wire [ADDR_WIDTH-1:0] addr_out
);

    localparam XOR_STAGES = TRNG_WIDTH / ADDR_WIDTH;
    
    wire [ADDR_WIDTH-1:0] stage_results [XOR_STAGES:0];
    
    assign stage_results[0] = addr_in;
    
    genvar i;
    generate
        for (i = 0; i < XOR_STAGES; i = i + 1) begin : xor_stages
            assign stage_results[i+1] = stage_results[i] ^ 
                   trng[(i+1)*ADDR_WIDTH-1 : i*ADDR_WIDTH];
        end
    endgenerate
    
    assign addr_out = stage_results[XOR_STAGES];

endmodule
