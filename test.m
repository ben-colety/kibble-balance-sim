clc
clear all

a = [51 72 33 44 55];
b = [21 82 43 34 25];
[a1, ai] = sort(a);
b1 = b(ai);

% a(2) = [];
% length(a)


%% structs
% syms alpha gamma
% p1 = struct('location', 'gnd-a1',   ...
%                 'type',     'cross',    ...
%                 'k',        4,          ...
%                 'ang_adm',  5,          ...
%                 'ener_var', alpha       );      %variable linked to energy
% p2 = struct('location', 'a1-a2',    ...
%                 'type',     'cross',    ...
%                 'k',        4,          ...
%                 'ang_adm',  5,          ...
%                 'ener_var', gamma       );
%             
% 
% pivots = [p1, p2];
% for i = 1 : length(pivots)
%     energies(i) = calcEnergy(pivots(i));
% end
%     
% function energy = calcEnergy(p)
%     energy = 0.5*p.k*p.ener_var^2;
% end

%%
% alpha = 5;
% beta = 5;
% alu = struct('location', 'a','angle_adm', 8, 'displacement', alpha+beta);
% if alu.displacement > alu.angle_adm
%     disp('broken')
% else
%     disp('fine')
% end

%%
% [x, y, z] = init(1, 2, 3);
% [x, y, z] = reset();
% 
% 
% 
% function [x, y, z] = init(a, b, c)
%     x = a;
%     y = b;
%     z = c;
% end
% 
% function [x, y, z] = reset()
%     x = 0;
%     y = 0;
%     z = 0;
% end