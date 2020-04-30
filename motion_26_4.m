clear all; clc;
%CONCEPTION DES MECANISMES II
%BALANCE DE KIBBLE PROJET 2020
%GROUP 11 SOLUTION 1
%BY Benjamin Colety
%
%SIMULATION OF FLEXIBLE DOUBLE WATT LINKAGE FOR KIBBLE BALANCE
%   this program simulates the motion and energy of our group's design of a
%   Kibble Balance moving from the weighted rest position (z = 0) to its
%   extremities in z. The balance implements flexible mechanisms which are
%   assumed to be perfectly machined. Vertical symmetry is assumed for the
%   behavior of the 2 watt linkages
%
%UNITS ARE (m, N, rad)

%% NOTATION

%   alpha   : angle between horizontal bar and horizon
%   beta    : angle between vertical bar and vertical
%   x       : horizontal displacment of end effector from origin
%   z       : vertical displacement of end effector from origin

%% materials initialization
safety_factor = 2;
alu = struct('E', 69e9, 'o_adm', 110e6/safety_factor);

materials = [alu];

%% position initialization
syms alpha beta z
gamma = alpha + beta;


%% watt linkage dimensions
L1 = 200e-3;            % length of horizontal bar
L2 = 200e-3;            % length of vertical bar
L3 = 250e-3;            % vertical distance between watt linkages
L4 = 300e-3;            % vertical distance between center of watt linkages
                        % and end effector
checkDims(L1, L2, L3, L4);

x_shift = 1e-6;

%% compensation dimensions
    %chariot
    Lc1 = 1;        % natural length of spring - initial compressed length
    Lc2 = 1;        % length of lames for parallel leaf spring
    Lc3 = 1;        % length of linking lame

    char_disp = Lc3 - sqrt(Lc3^2 - z^2);
    
    %weight
    
%% Pivots/Springs
%watt linkages
    p1 = struct('location', 'gnd-a1',   ...
                'type',     'cross',    ...
                'k',        4,          ...
                'cor_adm',  5,          ...
                'ener_var', alpha       );      %variable linked to energy
    p2 = struct('location', 'a1-a2',    ...
                'type',     'cross',    ...
                'k',        4,          ...
                'cor_adm',  5,          ...
                'ener_var', gamma       );
    p3 = p2;    p3.location = 'a2-a3';
    p4 = p1;    p4.location = 'a3-gnd';
    p5 = p1;    p5.location = 'gnd-b1';
    p6 = p2;    p6.location = 'b1-b2';
    p7 = p2;    p7.location = 'b2-b3';
    p8 = p1;    p8.location = 'b3-gnd';
    p9 = struct('location', 'a2-S',     ...
                'type',     'cross',    ...
                'k',        4,          ...
                'cor_adm',  5,          ...
                'ener_var', beta        );
    p10 = p9;   p10.location = 'b2-s';
    
%chariots
    sc = struct('location', 'gnd-c1',       ...
                'type',     'comp_spring',  ...
                'k',        4,              ...
                'cor_adm',  1000,           ... %N/A
                'ener_var', char_disp       );
    pc1 = struct('location', 'gnd-c1',      ...
                'type',     'parallel',     ...
                'k',        4,              ...
                'cor_adm',  5,              ...
                'ener_var', char_disp       );
    pc2 = struct('location', 'c1-c2',       ...
                'type',     'col',          ...
                'k',        4,              ...
                'cor_adm',  5,              ...
                'ener_var', asin(z/Lc3)     );
    pc3 = pc2; pc3.location = 'c2-s';
    
%weight comp
    sg = struct('location', 'gnd-s',        ...
                'type',     'spring',       ...
                'k',        4,              ...
                'cor_adm',  1000,           ... %N/A
                'ener_var', z               );
            
%% simulation
alpha_pas = 0.001;
alpha(1) = 0;
z(1) = 0;
x(1) = 0;
beta(1) = 0;
i = 1;
while abs(x(i)) < 1e-6
    plot(alpha(i), z(i), 'rx');hold on
    i=i+1;
    alpha(i) = alpha(i-1)+ alpha_pas;
    [beta(i), x(i), z(i)] = motionSim(alpha(i), L1, L2);
end
alpha(i) = 0;
z(i) = 0;
x(i) = 0;
beta(i) = 0;

figure(1);


%% FUNCTIONS
function checkDims(L1, L2, L3, L4)
    margin_x = 30e-3;
    margin_z = 200e-3;
    max_x = 600e-3;
    max_z = 800e-3;
    x_dim = 2*L1+margin_x;
    z_dim = L2/2 + L3/2 + L4 + margin_z;
    
    if  x_dim > max_x
        error('too big in x')
    elseif  z_dim > max_z
        error('too big in z')
    elseif  L3 < L2
        error('move linkages apart more (L3)')
    elseif  L4 < L3/2 + L2/2 + 50e-3
        error('make end effector farther away')
    end
end
function [beta, EE_x, EE_z] = motionSim(alpha, L1, L2)
    syms x z
    pt_2 = [L1*cos(alpha)-L1 L1*sin(alpha)-L2/2];
    pt_4 = [L1 L2/2];
    
    circ_2 = (x - pt_2(1))^2 + (z - pt_2(2))^2 == L2^2;
    circ_4 = (x - pt_4(1))^2 + (z - pt_4(2))^2 == L1^2;
    
    sols = solve(circ_2, circ_4, [x z]);
    
    if  sols.z(1) > pt_2(2) && sols.z(1) > sols.z(2) 
        pt_3 = [sols.x(1) sols.z(1)];
    elseif sols.z(2) > pt_2(2) && sols.z(2) > sols.z(1) 
        pt_3 = [sols.x(2) sols.z(2)];
    else
        error('couldnt solve for pt_3');
    end
    
    beta = atan((pt_3(1)-pt_2(1))/(pt_3(2)-pt_2(2)));
    EE = pt_2 + L2/2*[sin(beta) cos(beta)];
    EE_x = EE(1);
    EE_z = EE(2);
end


