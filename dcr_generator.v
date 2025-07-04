`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.06.2025 14:46:19
// Design Name: 
// Module Name: dcr_generator
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


module dcr_generator (
    input wire clk,
    input wire rst_n,
    input wire boot_event,
    input wire power_event,
    input wire tamper_event,
    output wire dcr
);

    assign dcr = boot_event | power_event | tamper_event | ~rst_n;

endmodule
