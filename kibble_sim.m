
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
cast_alu = struct('E', 72e9, 'o_adm', 200e6/safety_factor);
wrought_alu = struct('E', 72e9, 'o_adm', 330e6/safety_factor);
brass = struct('E', 100e9, 'o_adm', 200e6/safety_factor);
bronze = struct('E', 100e9, 'o_adm', 400e6/safety_factor);
cast_mag = struct('E', 45e9, 'o_adm', 150e6/safety_factor);
titan = struct('E', 100e9, 'o_adm', 300e6/safety_factor);

Materials = [cast_alu wrought_alu brass bronze cast_mag titan];

%% Position Initialization
syms alpha beta z phi
gamma = alpha + beta;

%% Physical Dimensions

    %General
    m_pesee = 3;
    m_system = 3;
    m = m_pesee + m_system;
    g = 9.81; 
    b = 50e-3;              %thickness of planar material

    %Linear Acuator dimensions
    La1 = 114.3e-3;         %motor length at mid-stroke
    La2 = 31.8e-3;          %diameter of motor housing
    
    
    %Watt Linkage
    L1 = 250e-3;            % length of horizontal bar
    L2 = 75e-3;             % length of vertical bar
    L3 = 250e-3;            % vertical distance between watt linkages
    L4 = 300e-3;            % vertical distance between center of "S"
                            % and end effector

    %Rigidity Compensation
    Lc1 = 50e-3;            % distance spring is compressed when z = 0
    Lc2 = 30e-3;            % length of lames for parallel leaf spring
    Lc3 = 100e-3;           % length of linking lame

    char_disp = Lc3 - sqrt(Lc3^2 - z^2);
    
    %Weight Compensation
    Lg = m*g/1000;          %extension of spring at z = 0 for k = 1000
    
    checkDims(L1, L2, L3, L4, La1, La2, Lc1, Lc2, Lc3, Lg, b);
    
%% Pivots/Springs
%watt linkages
    p1 = pivot('point', 4, alpha, 4);          %type, k, displacement variable
    p2 = pivot('point', 3, gamma, 4);  
    p3 = pivot('point', 1, beta,  2);
    
    
%rigidity comp
    sc  = pivot('spring', 80, Lc1-char_disp, 2, 3);
    pc1 = pivot('parallel', 5, char_disp, 2);
    pc2 = pivot('point', 5, asin(z/Lc3), 4);
    
    
%weight comp
    sg  = pivot('spring', 1, z+Lg, 1, 3    );
            
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

gamma = eval(gamma);
Positions = [alpha beta gamma z x];

%% Pivot/Spring Physical Dimensions

%for modifying k's easily
p1.k = p1.k;
p2.k = p2.k;
p3.k = p3.k;
sc.k = sc.k;
pc1.k = pc1.k;
pc2.k = pc2.k;
sg.k = sg.k;
        
for i = 1:length(Materials)            
    Pivots(:,i) = [p1, p2, p3, sc, pc1, pc2, sg];
end

for i = 1:length(Materials)
    for j = 1:length(Pivots(:,i))
        Pivots(j,i).verified = [1 1];
        switch Pivots(i).type
            case {'spring','parallel'}
                syms h L
                rig = Pivots(j,i).k == Pivots(j,i).num_lames * Materials(i).E * b * h^3 / L^3;
                adm = max(abs(eval(Pivots(j,i).ener_var))) == Materials(i).o_adm*L^2 /(3*Materials(i).E*h);
                Pivots(j,i) = pivotSolve(Pivots(j,i),rig, adm, h, L);
                %dimension verification
                if Pivots(j,i).h < 100e-6
                    warning('Pivot %d of type %s in material %d has blades that are too thin\n',j,Pivots(j,i).type,i)
                    Pivots(j,i).verified(1) = 0;
                elseif Pivots(j,i).verified(1) == 1
                    fprintf('Pivot %d of type %s in material %d seems to work but the blade length still needs to be checked\n',j,Pivots(j,i).type,i)
                end

            case {'point','col','cross'}
                %calculation as a col
                syms e r
                rig = Pivots(j,i).k == 2* Materials(i).E * b * e^(2.5) / (9*pi*r^(0.5));
                adm = max(abs(eval(Pivots(j,i).ener_var))) == 3*pi*Materials(i).o_adm*sqrt(r)/(4*Materials(i).E*sqrt(e));
                Pivots(j,i) = pivotSolve(Pivots(j,i),rig, adm, e, r);
                %dimension verification
                if Pivots(j,i).r/Pivots(j,i).e < 5
                    warning('Pivot %d as a col in Material %d \tr/e ratio too large', j, i)
                    Pivots(j,i).verified(1) = 0;
                end
                if Pivots(j,i).e < 100e-6 % moodle -> b < 20mm + de passage de fil -> e ~=50um
                    warning('Pivot %d as a col in Material %d \te too small : %d', j, i, Pivots(j,i).e)
                    Pivots(j,i).verified(1) = 0;
                elseif Pivots(j,i).e > 1e-3
                    warning('Pivot %d as a col in Material %d \te too large : %d', j, i, Pivots(j,i).e)
                    Pivots(j,i).verified(1) = 0;
                end
                if Pivots(j,i).r < 1e-4
                    warning('Pivot %d as a col in Material %d \tr too small : %d', j, i, Pivots(j,i).r)
                    Pivots(j,i).verified(1) = 0;
                elseif Pivots(j,i).r > 0.01 %r must be smaller than 1cm
                    warning('Pivot %d as a col in Material %d \tr too big : %d', j, i, Pivots(j,i).r)
                    Pivots(j,i).verified(1) = 0;
                end
                if Pivots(j,i).verified(1) == 1
                    fprintf('Pivot %d as a col in Material %d with e : %d, r : %d is acceptable \n', j, i, Pivots(j,i).e, Pivots(j,i).r);
                end

                %calculation as a cross
                syms h L
                rig = Pivots(j,i).k == 8*Materials(i).E*b*h^3 /(12*L);
                adm = max(abs(eval(Pivots(j,i).ener_var))) == Materials(i).o_adm * L /(2*Materials(i).E*h);
                Pivots(j,i) = pivotSolve(Pivots(j,i), rig, adm, h, L);
                %dimension verification
                if Pivots(j,i).L/Pivots(j,i).h > 60
                    warning('Pivot %d as a cross in Material %d \tL/h ratio too large', j, i)
                    Pivots(j,i).verified(2) = 0;
                end
                if Pivots(j,i).h < 100e-6
                    warning('Pivot %d as a cross in Material %d \th too large : %d', j, i, Pivots(j,i).h)
                    Pivots(j,i).verified(2) = 0;
                end
                if Pivots(j,i).verified(2) == 1
                    fprintf('Pivot %d as a cross in Material %d has good dimensions \n', j, i)
                end
        end
    end
end

%internal forces check
    %rigidity compensation
    %pc2 as col
    for i = 1:length(Materials)
        spring_force = 2*Pivots(4,i).k*Lc1;
        area = Pivots(6).e*b;
        fail = (spring_force/area) >= Materials(i).o_adm;
        if fail
            fprintf('Pivot Pc1 can''t take the pressure with material %d',i);
        end
    end

%% Energy & Force
Energies = zeros(length(z), length(Pivots)+1);
Energies(:,length(Pivots)+1) = -m*g*z;
    for i = 1:length(Pivots)
        ener = Pivots(i,1).quant*calcEnergy(Pivots(i,1));
        Energies(:,i) = eval(ener);
    end

sumEnergies = sum(Energies, 2);     %summing energies of pivots for each z

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
yline(1e-6, 'k', '+1\mum');
yline(-1e-6, 'k', '-1\mum');

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
function checkDims(L1, L2, L3, L4, La1, La2, Lc1, Lc2, Lc3, Lg, b)
    margin_x = 30e-3;
    margin_z = 30e-3;
    max_x = 2*sqrt(((600e-3)/2)^2 - (b/2)^2);       %taking into account cylindrical constraint
    max_z = 800e-3;
    x_dim = 2*L1+2*margin_x;
    z_dim = L2/2 + L3/2 + L4 + 2*margin_z;
    
    if  x_dim > max_x
        error('too big in x')
    elseif  z_dim > max_z
        error('too big in z')
    elseif  L3 < L2
        error('move linkages apart more (L3)')
    elseif  L4 < L3/2 + L2/2 + margin_z             %can also add La1 if motor is mounted below linkages
        error('make end effector farther away')
    elseif  L3-L2 < La1 + margin_z
        error('not enough space to mount motor inside linkages: %.2f > %.2f', La1 + margin_z, L3-L2)
    elseif  max_x < 2*Lc3 + 2*30e-3 + 2*2*Lc1+2*margin_x
        error('rigidity compensation too wide: %.2f > %.2f',2*Lc3 + 2*30e-3 + 2*2*Lc1+2*margin_x,max_x)
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

