
%CONCEPTION DES MECANISMES II
%BALANCE DE KIBBLE PROJET 2020
%GROUP 11 SOLUTION 1
%BY Benjamin Colety
%
%SIMULATION OF WATT LINKAGE FOR KIBBLE BALANCE
%   this program simulates the motion of our group's design of a
%   Kibble Balance moving from the weighted rest position (z = 0) to its
%   extremities in z. The balance implements flexible mechanisms which are
%   assumed to be perfectly machined. Vertical symmetry is assumed for the
%   behavior of the 2 watt linkages
%
%UNITS ARE (m, N, rad)

%% NOTATION

%   alpha   : angle between horizontal bar and horizon
%   beta    : angle between vertical bar and vertical
%   gamma1  : angle between horizontal bar and vertical bar(1)
%   gamma2  : angle between horizontal bar and vertical bar(2)
%   x       : horizontal displacment of end effector from origin
%   z       : vertical displacement of end effector from origin
clear all; clc; close all; format shorteng;

%% Position Initialization
syms alpha beta z
gamma1 = alpha + beta;
gamma2 = alpha - beta;

%% Physical Dimensions

    %General
    m_pesee = 3;
    m_system = 3;
    m = m_pesee + m_system;
    g = 9.81; 
    b = 30e-3;              %thickness of planar material
    
    %Watt Linkage
    L1 = 250e-3;            % length of horizontal bar
    L2 = 75e-3;             % length of vertical bar

%% Simulation

%initialization
alpha_pas = 0.001;
alpha = zeros(500, 1);
beta = zeros(500, 1);
x = zeros(500, 1);
z = zeros(500, 1); 
i = 1;
%linkage descending
while (abs(x(i)) < 1e-6 && abs(z(i)) < 15e-3) || abs(z(i)) < 16e-3
    i=i+1;
    alpha(i) = alpha(i-1)+ alpha_pas;
    [beta(i), x(i), z(i)] = motionSim(alpha(i), L1, L2);
end
%reset to 0
i=i+1;
%rising
while (abs(x(i)) < 1e-6 && abs(z(i)) < 15e-3) || abs(z(i)) < 16e-3
    i=i+1;
    alpha(i) = alpha(i-1)- alpha_pas;
    [beta(i), x(i), z(i)] = motionSim(alpha(i), L1, L2);
end

%trimming arrays
alpha = alpha(2:find(alpha,1,'last'));
beta = beta(2:find(beta,1,'last'));
x = x(2:find(x,1,'last'));
z = z(2:find(z,1,'last'));

%sorting arrays
[alpha, sort_ind] = sort(alpha);
beta = beta(sort_ind);
z = z(sort_ind);
x = x(sort_ind);

gamma1 = eval(gamma1);
gamma2 = eval(gamma2);
Positions = [alpha beta gamma1 gamma2 z x];

%% Results

max_alpha   = max(abs(alpha))
max_beta    = max(abs(beta))
max_gamma   = max(abs(gamma1))

disp('during the linear movement');
z_course_ind = find(abs(z) < 15e-3);
%maximum parasitic motion in x during linear trajectory
max_x       = max(abs(x(min(z_course_ind):max(z_course_ind))))

%% Graphics
figure(1);
plot1 = plot(z, x, 'r');
title('z vs x');
xlabel('z (m)');
ylabel('x (m)');
xline(15e-3,'k','+15mm');
xline(-15e-3,'k','-15mm');
yline(1e-6, 'k', '+1\mum');
yline(-1e-6, 'k', '-1\mum');


%% FUNCTIONS
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

