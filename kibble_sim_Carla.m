
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
clear all; clc;

%% materials initialization
safety_factor = 2;
alu = struct('E', 69e9, 'o_adm', 110e6/safety_factor);

Materials = [alu];

%% Position Initialization
syms alpha beta z phi
gamma = alpha + beta;

%% Physical Dimensions
    %Watt Linkage
    L1 = 150e-3;            % length of horizontal bar
    L2 = 75e-3;             % length of vertical bar
    L3 = 250e-3;            % vertical distance between watt linkages
    L4 = 300e-3;            % vertical distance between center of "S"
                            % and end effector

    %Mass
    m = 3;                              %Need to add Linkage name
    g = 9.81;
    b = 30e-3;              % system width
    
    %Rigidity Compensation
    Lc1 = 50e-3;            % distance spring is compressed when z = 0
    Lc2 = 30e-3;            % length of lames for parallel leaf spring
    Lc3 = 50e-3;            % length of linking lame

    char_disp = Lc3 - sqrt(Lc3^2 - z^2);
    
    %Weight Compensation
    Lg = m*g/1000;          %extension of spring at z = 0 for k = 1000
    
    %Lever
    Llev1 = 150e-3;               %full length of lever
    Llev2 = 85e-3;               %length of section from pivot to application point
    phi = asin(z/Llev2);    %angle between lever and horizontal
    zm = Llev1*sin(phi);    %position of motor
    
    
    
    checkDims(L1, L2, L3, L4);
    
%% Pivots/Springs
%watt linkages
    p1 = struct('location', 'gnd-a1',   ...
                'type',     'cross',    ...
                'k',        1,        ...     %effective rigidity
                'cor_adm',  5,          ...     %course admissible
                'ener_var', alpha,      ...     %variable linked to energy
                'Dims',     zeros(1,8)  );      %L, h, e/r, e_min, e_max, r_min, r_max, d_compression/ pivots croix et table lame parallel ; pivots cols ; ressorts
    p2 = struct('location', 'a1-a2',    ...
                'type',     'cross',    ...
                'k',        1,        ...
                'cor_adm',  5,          ...
                'ener_var', gamma,      ...
                'Dims',     zeros(1,8)  );
    p3 = p2;    p3.location = 'a2-a3';
    p4 = p1;    p4.location = 'a3-gnd';
    p5 = p1;    p5.location = 'gnd-b1';
    p6 = p2;    p6.location = 'b1-b2';
    p7 = p2;    p7.location = 'b2-b3';
    p8 = p1;    p8.location = 'b3-gnd';
    p9 = struct('location', 'a2-S',     ...
                'type',     'cross',    ...
                'k',        1,       ...
                'cor_adm',  5,          ...
                'ener_var', beta,       ...
                'Dims',     zeros(1,8)  );
    p10 = p9;   p10.location = 'b2-s';
    
%rigidity comp
    sc = struct('location',  'gnd-c1',      ...
                 'type',     'comp_spring', ...
                 'k',        1,          ...
                 'cor_adm',  1000,          ... %N/A
                 'ener_var', Lc1-char_disp, ...
                 'Dims',     zeros(1,8)  );
    pc1 = struct('location', 'gnd-c1',      ...
                 'type',     'parallel',    ...
                 'k',        1,           ...
                 'cor_adm',  5,             ...
                 'ener_var', char_disp,     ...
                 'Dims',     zeros(1,8)     );
    pc2 = struct('location', 'c1-c2',       ...
                 'type',     'col',         ...
                 'k',        1,          ...
                 'cor_adm',  5,             ...
                 'ener_var', asin(z/Lc3),   ...
                 'Dims',     zeros(1,8)     );
    pc3 = pc2; pc3.location = 'c2-s';
    
%weight comp
    sg = struct('location', 'gnd-s',        ...
                 'type',     'spring',      ...
                 'k',        1,          ...
                 'cor_adm',  1000,          ... %N/A
                 'ener_var', z+Lg,          ...
                 'Dims',     zeros(1,8)     );
            
%lever
    pl1 = struct('location', 'gnd-lev2',    ...
                 'type',     'parallel',    ...
                 'k',        1,          ...
                 'cor_adm',  5,             ...
                 'ener_var', zm,            ...
                 'Dims',     zeros(1,8)     );             %need to solve
    pl2 = struct('location', 'lev2-lev1',   ...
                 'type',     'lame',        ...
                 'k',        1,          ...
                 'cor_adm',  5,             ...
                 'ener_var', phi,           ...
                 'Dims',     zeros(1,8)     );             %need to solve
    pl3 = struct('location', 'lev1-s',      ...
                 'type',     'lame',        ...
                 'k',        1,          ...
                 'cor_adm',  5,             ...
                 'ener_var', phi,           ...
                 'Dims',     zeros(1,8)     );              %need to solve
    pl4 = struct('location', 'lev1-gnd',    ...
                 'type',     'col',         ...
                 'k',        1,          ...
                 'cor_adm',  5,             ...
                 'ener_var', phi,           ...
                 'Dims',     zeros(1,8)     );              %need to solve
            
Pivots_link = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10];
Pivots_rigidity = [sc, pc1, pc2, pc3];
Pivots_weight = [sg];
Pivots_lever = [pl1, pl2, pl3, pl4];
Pivots = [Pivots_link, Pivots_rigidity, Pivots_weight, Pivots_lever];
            
%% Simulation

%initialization
% alpha_pas = 0.001;
% alpha = zeros(500, 1);
% beta = zeros(500, 1);
% x = zeros(500, 1);
% z = zeros(500, 1);
% beta = zeros(500, 1);
% i = 1;
% %linkage descending
% while (abs(x(i)) < 1e-6 && abs(z(i)) < 15e-3) || abs(z(i)) < 30e-3
%     i=i+1;
%     alpha(i) = alpha(i-1)+ alpha_pas;
%     [beta(i), x(i), z(i)] = motionSim(alpha(i), L1, L2);
% end
% %reset to 0
% i=i+1;
% %rising
% while (abs(x(i)) < 1e-6 && abs(z(i)) < 15e-3) || abs(z(i)) < 30e-3
%     i=i+1;
%     alpha(i) = alpha(i-1)- alpha_pas;
%     [beta(i), x(i), z(i)] = motionSim(alpha(i), L1, L2);
% end
% 
% %trimming arrays
% alpha = alpha(2:find(alpha,1,'last'));
% beta = beta(2:find(beta,1,'last'));
% x = x(2:find(x,1,'last'));
% z = z(2:find(z,1,'last'));
% 
% %sorting arrays
% [alpha, sort_ind] = sort(alpha);
% beta = beta(sort_ind);
% z = z(sort_ind);
% x = x(sort_ind);
% max(abs(alpha))
% max(abs(beta))
% %max(abs(gamma))
% %max(abs(phi))
% max(abs(z))
% max(abs(x))
alpha_m=0.2030; beta_m=0.0815; z_m=0.0301; x_m=2.5539e-05;

%% Energy & Force
% Energies = zeros(length(z), length(Pivots));
% for i = 1:length(Pivots)
%     ener = calcEnergy(Pivots(i));
%     Energies(:,i) = eval(ener);
% end
% sumEnergies = sum(Energies, 2);
% 
% %Force
% sampling_res = max(z)*2/(length(z)*4);
% sampling_pts = min(z):sampling_res:max(z);
% fctEnergy = interp1(z,sumEnergies, sampling_pts, 'spline');
% force = diff(fctEnergy);
% 
% %% Graphics
% figure(1);
% plot(z, x, 'r');
% 
% figure(2);
% plot(z, sumEnergies, 'r'); hold on %plot sum of all Energies
% figure(3);
% plot(sampling_pts(1:length(sampling_pts)-1), force, 'g'); %plot Forces

%% Dimensioning Pivots
for i=1:length(Pivots)
    j=1;                                % choix matériaux = alu
    if ((Pivots(i).type == "cross") | (Pivots(i).type == "col")) % check if pivots
        h_L = Materials(j).o_adm/(2*(Materials(j).E)*alpha_m); %change recording to max of angle used
        hmax = sqrt(3*Pivots(i).k/(2*(Materials(j).E)*b*h_L));
        Lmin = hmax/h_L;
        Pivots(i).Dims(1) = Lmin; Pivots(i).Dims(2) = hmax; % add in table L h for cross
        
%         syms e r
%         f = Pivots(i).k == 2*Materials(j).E*b*e^2.5/pi/sqrt(r);
%         g = alpha_m == 3*pi*Materials(j).o_adm*sqrt(r)/4/Materials(j).E/sqrt(e);
%         solutions = solve(f,g,[e,r])
        e_r = (3*pi*Materials(j).o_adm/(4*Materials(j).E*alpha_m))^2;
        if e_r > 1/5
            warning('Pivot à col %d has a too large e_r ratio', i)
        end
        e_max= b/10; 
        if e_max > 1e-3
            e_max = 1e-3;
            e_min = 1e-6;
        elseif e_max < 1e-6
            e_max = -1e-6; e_min = -1e-6;
            warning('Pivot à col %d has a too large e_min', i)
        else
            e_min = 1e-6;
        end
        %L, h, e/r, e_min, e_max, r_min, r_max, d_compression
        r_min = e_max/e_r;
        r_max = e_min/e_r;
%         if r_min > r_max
%             warning('Pivot à col %d has r_min > r_max', i)
%             r_min = -1e-6; r_max = -1e-6;
%         end
        if r_max > 1
            r_max = 1;
        elseif r_max < 1e-4
            warning('Pivot à col %d has too small r_max', i)
            r_min = -1e-6; r_max = -1e-6;
        end

        if r_min < 1e-4
            r_min = 1e-4;
        elseif r_min > 1
            warning('Pivot à col %d has r_min too big', i)
        end
        Pivots(i).Dims(3) = e_r; Pivots(i).Dims(4) = e_max; 
        Pivots(i).Dims(5) = e_min; Pivots(i).Dims(6) = r_min; 
        Pivots(i).Dims(3) = r_max;
        end  
        Pivots(i).Dims
%         
   % elseif 
        
%     elseif p.type == 'comp_spring'
%         Lmin = 'not defined';
%         hmax = Lmin;
%     elseif p.type == 'parallel'
%         Lmin = 'not defined';
%         hmax = Lmin;
%     elseif p.type == 'spring'
%         Lmin = 'not defined';
%         hmax = Lmin;
    end

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

