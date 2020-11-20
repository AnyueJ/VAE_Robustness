% This function is written to calculate the pressure field for various 
% Pumping tests based on the location of Pumping test
% The output is a matrix with collumns coressponding to the pressures of
% each index

function [P_OUT] = pressure_calculation_main(Hyd_Con);


P_OUT = [];

    

    p = pressure_calculation(Hyd_Con);
    P_OUT = [P_OUT p];

