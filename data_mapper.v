`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.06.2025 14:46:19
// Design Name: 
// Module Name: data_mapper
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


module data_mapper #(
    parameter DATA_WIDTH = 52,
    parameter TRNG_WIDTH = 32
)(
    input wire [DATA_WIDTH-1:0] data_in,
    input wire [TRNG_WIDTH-1:0] trng,
    output wire [DATA_WIDTH-1:0] data_out
);

    // Calculate number of complete XOR stages
    localparam XOR_STAGES = (DATA_WIDTH + TRNG_WIDTH - 1) / TRNG_WIDTH;
    
    wire [DATA_WIDTH-1:0] stage_results [XOR_STAGES:0];
    
    assign stage_results[0] = data_in;
    
    genvar i;
    generate
        for (i = 0; i < XOR_STAGES; i = i + 1) begin : xor_stages
            if ((i+1) * TRNG_WIDTH <= DATA_WIDTH) begin
                // Full TRNG width XOR
                assign stage_results[i+1][(i+1)*TRNG_WIDTH-1 : i*TRNG_WIDTH] = 
                       stage_results[i][(i+1)*TRNG_WIDTH-1 : i*TRNG_WIDTH] ^ trng;
                
                // Pass through remaining bits unchanged
                if (i*TRNG_WIDTH > 0) begin
                    assign stage_results[i+1][i*TRNG_WIDTH-1 : 0] = 
                           stage_results[i][i*TRNG_WIDTH-1 : 0];
                end
                if ((i+1)*TRNG_WIDTH < DATA_WIDTH) begin
                    assign stage_results[i+1][DATA_WIDTH-1 : (i+1)*TRNG_WIDTH] = 
                           stage_results[i][DATA_WIDTH-1 : (i+1)*TRNG_WIDTH];
                end
            end else begin
                // Partial TRNG width XOR for last stage
                localparam REMAINING_BITS = DATA_WIDTH - i*TRNG_WIDTH;
                assign stage_results[i+1][DATA_WIDTH-1 : i*TRNG_WIDTH] = 
                       stage_results[i][DATA_WIDTH-1 : i*TRNG_WIDTH] ^ 
                       trng[REMAINING_BITS-1 : 0];
                
                // Pass through lower bits unchanged
                if (i*TRNG_WIDTH > 0) begin
                    assign stage_results[i+1][i*TRNG_WIDTH-1 : 0] = 
                           stage_results[i][i*TRNG_WIDTH-1 : 0];
                end
            end
        end
    endgenerate
    
    assign data_out = stage_results[XOR_STAGES];

endmodule
