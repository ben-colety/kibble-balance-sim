clc; clear all;

%Material
%wrought alu
E = 72e9;
o_adm = 165e6;

%Thickness
b = 30e-3;
%Lames
num_lames = ;
h = ;
L = ;

rig = num_lames * E * b * h^3 / L^3;
adm = o_adm * L^2 /(3 * E * h);

