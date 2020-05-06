clc 
clear all

syms x y

f = (x*y) == 1;
g = (x)+y == 0;
huh = [{x} {y}];
solution = testfunct(f, g, {x}, {y})


function sols1 = testfunct(eq1, eq2, t, s)

    sols = solve(eq1, eq2, [t s]);
    sols1 = sols.((t{1,1}))
%     sols1 = sols.t(1)(imag(sols.t)==0);
%     sols2 = sols.s(1)(imag(sols.s)==0);
end
%%
% phi = 0:.01:5;
% Positions = phi;
% alu = struct('E', 69e9, 'o_adm', 110e6/2);
% 
% Materials = [alu];
% b = 50e-3;
%  pl2 = struct('location', 'lev2-lev1',   ...
%                  'type',     'lame',        ...
%                  'k',        1,             ...
%                  'cor_adm',  5,             ...
%                  'ener_var', phi,           ... 
%                  'Dims',     [3; zeros(7,1)]);
%  pl4 = struct(   'location', 'lev1-gnd',    ...
%                  'type',     'col',         ...
%                  'k',        0.5,           ...
%                  'cor_adm',  5,             ...
%                  'ener_var', phi,           ... 
%                  'Dims',     zeros(1, 8)    );
% 
% Pivots = [pl2 pl4];
% 
% syms h L
% for i = 1:2;
% for j = 1;
% switch Pivots(i).type
%     case 'lame'
%         rig = Pivots(i).k == Pivots(i).Dims(1) * Materials(j).E * b * h^3 / L^3;
%         adm = max(abs(Pivots(i).ener_var)) == Materials(j).o_adm*L^2 /(3*Materials(j).E*h);
% 
% %         solutions = solve(rig, adm, [h L]);
% % 
% %         h_sol  = eval(solutions.h)
% %         L_sol  = eval(solutions.L)
% %         h_real = eval(solutions.h(imag(solutions.h)==0))
% %         L_real = eval(solutions.L(imag(solutions.L)==0))
%         possible = pivotSolve(rig, adm, {h},{L})
% 
% end
% end
% end
% 
% function solutions = pivotSolve(eq1, eq2, var1, var2)
%         
%         solutions = solve(eq1, eq2);
% 
%         h_sol  = eval(solutions.var1)
%         L_sol  = eval(solutions.var2)
%         h_real = eval(solutions.var1(imag(solutions.var1)==0))
%         L_real = eval(solutions.var2(imag(solutions.var2)==0))
% 
%         pos_solutions = struct('h', h_real, 'L', L_real);
% end

%%
% a = struct( 'x', 1, 'y', 2, 'dims', zeros(1,8));
% b = struct( 'x', 1, 'y', 2, 'dims', zeros(1,8));
% 
% tab = [a b];
% 
% tab = setx(tab, 1, 5);
% tab(1)
% 
% function  tab = setx(tab, i, t)
%     tab(i).x = t;
%     tab(i).dims(5) = t;
% end

%%
% syms x y
%    sc = struct(  'location',  'gnd-c1',      ...
%                  'type',     'comp_spring', ...
%                  'k',        80,            ...
%                  'cor_adm',  1000,          ... %N/A
%                  'ener_var', x  );
%     pc1 = struct('location', 'gnd-c1',      ...
%                  'type',     'parallel',    ...
%                  'k',        5,             ...
%                  'cor_adm',  5,             ...
%                  'ener_var', x      );
%     pc2 = struct('location', 'c1-c2',       ...
%                  'type',     'col',         ...
%                  'k',        5,             ...
%                  'cor_adm',  5,             ...
%                  'ener_var', y    );
% 
%              
%      pivots = [sc pc1 pc2];
%      T = struct2table(pivots)
%      T(4) = 

%%
% a = "help"
% b = "help"
% c = "no help"
% 
% if (a == b) || (b == c)
%     disp('help is on the way');
% end

%%

% a = [51 72 33 44 55];
% b = [21 82 43 34 25];
% [a1, ai] = sort(a);
% b1 = b(ai);

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