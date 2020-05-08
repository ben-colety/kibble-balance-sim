
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
clear all; clc; close all; format shorteng;

%% materials initialization
safety_factor = 2;
alu = struct('E', 69e9, 'o_adm', 110e6/safety_factor);
steel = struct('E', 45, 'o_adm', 120e6/safety_factor);

Materials = [alu steel];

%% Position Initialization
syms alpha beta z phi
gamma = alpha + beta;

%% Physical Dimensions

    %General
    m = 6;                          %Need to add Linkage mass
    g = 9.81;
    b = 50e-3;              %thickness of planar material

    %Watt Linkage
    L1 = 250e-3;            % length of horizontal bar
    L2 = 75e-3;             % length of vertical bar
    L3 = 250e-3;            % vertical distance between watt linkages
    L4 = 300e-3;            % vertical distance between center of "S"
                            % and end effector

    %Rigidity Compensation
    Lc1 = 50e-3;            % distance spring is compressed when z = 0
    Lc2 = 30e-3;            % length of lames for parallel leaf spring
    Lc3 = 200e-3;           % length of linking lame

    char_disp = Lc3 - sqrt(Lc3^2 - z^2);
    
    %Weight Compensation
    Lg = m*g/1000;          %extension of spring at z = 0 for k = 1000
    
    checkDims(L1, L2, L3, L4);
    
%% Pivots/Springs
%watt linkages
    p1 = pivot('point', 4, alpha, 4);          %type, k, displacement variable
    p2 = pivot('point', 3, gamma, 4);  
    p3 = pivot('point', 1, beta,  2);
    
    Pivots_link = [p1, p2, p3];
    
%rigidity comp
    sc  = pivot('spring', 80, Lc1-char_disp, 2, 3);
    pc1 = pivot('parallel', 5, char_disp, 2);
    pc2 = pivot('point', 5, asin(z/Lc3), 4);
    
    Pivots_rigidity = [sc, pc1, pc2];
    
%weight comp
    sg  = pivot('spring', 1, z+Lg, 1, 3    );
             
    Pivots_weight = [sg];
for i = 1:length(Materials)            
Pivots(:,i) = [Pivots_link, Pivots_rigidity, Pivots_weight];
end
            
%% Simulation

%initialization
alpha_pas = 0.001;
alpha = zeros(500, 1);
beta = zeros(500, 1);
x = zeros(500, 1);
z = zeros(500, 1); 
i = 1;
%linkage descending
while (abs(x(i)) < 1e-6 && abs(z(i)) < 15e-3) || abs(z(i)) < 30e-3
    i=i+1;
    alpha(i) = alpha(i-1)+ alpha_pas;
    [beta(i), x(i), z(i)] = motionSim(alpha(i), L1, L2);
end
%reset to 0
i=i+1;
%rising
while (abs(x(i)) < 1e-6 && abs(z(i)) < 15e-3) || abs(z(i)) < 30e-3
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

gamma = eval(gamma);
Positions = [alpha beta gamma z x];

%% Pivot/Spring Physical Dimensions

for i = 1:length(Materials)
    for j = 1:length(Pivots(:,i))
        switch Pivots(i).type
            case {'spring','parallel'}
                syms h L
                rig = Pivots(j,i).k == Pivots(j,i).num_lames * Materials(i).E * b * h^3 / L^3;
                adm = max(abs(eval(Pivots(j,i).ener_var))) == Materials(i).o_adm*L^2 /(3*Materials(i).E*h);
                if Pivots(j,i).type == "parallel"
                    fprintf('parallel for pivot %d\n',j)
                elseif Pivots(j,i).type == "spring"
                    fprintf('spring for pivot %d with %d lames\n',j,Pivots(j,i).num_lames)
                end
                Pivots(j,i) = pivotSolve(Pivots(j,i),rig, adm, h, L)

            case {'point','col','cross'}
                %calculation as a col
                syms e r
                rig = Pivots(j,i).k == 2* Materials(i).E * b * e^(2.5) / (9*pi*r^(0.5));
                adm = max(abs(eval(Pivots(j,i).ener_var))) == 3*pi*Materials(i).o_adm*sqrt(r)/(4*Materials(i).E*sqrt(e));
                fprintf('col for pivot %d\n',j)
                Pivots(j,i) = pivotSolve(Pivots(j,i),rig, adm, e, r)

                %calculation as a cross
                syms h L
                rig2 = Pivots(j,i).k == 8*Materials(i).E*b*h^3 /(12*L);
                adm2 = max(abs(eval(Pivots(j,i).ener_var))) == Materials(i).o_adm * L /(2*Materials(i).E*h);
                fprintf('cross for pivot %d\n',j)
                Pivots(j,i) = pivotSolve(Pivots(j,i), rig2, adm2, h, L)
        end
    end
end

%internal forces check



%% Energy & Force
Energies = zeros(length(z), length(Pivots)+1);
for i = 1:length(Pivots)
    ener = Pivots(i).quant*calcEnergy(Pivots(i));
    Energies(:,i) = eval(ener);
end
Energies(:,length(Pivots)+1) = -m*g*z;
sumEnergies = sum(Energies, 2);

%Force
sampling_res = max(z)*2/(length(z)*4);
sampling_pts = min(z):sampling_res:max(z);
fctEnergy = interp1(z,sumEnergies, sampling_pts, 'spline');
force = diff(fctEnergy);

%Rigidity Tangentielle Residuelle
rig_res = diff(force);



%% Results

max_alpha   = max(abs(alpha))
max_beta    = max(abs(beta))
max_gamma   = max(abs(gamma))

disp('during the linear movement');
z_course_ind = find(abs(z) < 15e-3);
tmp = max(z_course_ind);
%maximum parasitic motion in x during linear trajectory
max_x       = max(abs(x(min(z_course_ind):max(z_course_ind))))
z_course_ind = find(abs(sampling_pts) < 15e-3);
%maximum residual force during linear trajectory
max_force   = max(abs(force(min(z_course_ind):max(z_course_ind))))


%% Graphics
%important results
max_z_graph = tmp;


figure(1);
plot1 = plot(z, x, 'r');
title('z vs x');
xlabel('z (m)');
ylabel('x (m)');
xline(15e-3,'k','+15mm');
xline(-15e-3,'k','-15mm');

figure(2);
plot(z, sumEnergies, 'r'); %plot sum of all Potential Energies
title('z vs Sum of All Potentiel Energies')
xlabel('z (m)');
ylabel('Energy (J)');
xline(15e-3,'k','+15mm');
xline(-15e-3,'k','-15mm');

figure(3);
plot(sampling_pts(1:length(sampling_pts)-1), force, 'g'); %plot Forces
title('z vs Residual Force')
xlabel('z (m)');
ylabel('Residual Force (N)');
xline(15e-3,'k','+15mm');
xline(-15e-3,'k','-15mm');

figure(4);
plot(sampling_pts(1:length(sampling_pts)-2),rig_res);
title('z vs Residual Tangential Rigidity')
xlabel('z (m)');
ylabel('Residual Rigidity (N/m)');
xline(15e-3,'k','+15mm');
xline(-15e-3,'k','-15mm');

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
function energy = calcEnergy(p)
    energy = 0.5*p.k*p.ener_var^2;
end
function stop = checkTolerances(Materials, Pivots, alpha, beta, z, x)
    %parasitic motion tolerances
    limit_trans = 1e-6;
    limit_rot_xy = 50e-6;
    limit_rot_z = 1e-6;
    
end
function failure = checkFlambage(Pivots, Materials, Positions)

end

