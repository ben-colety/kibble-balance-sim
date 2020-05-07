clc 
clear all

% syms g h
% b_1 = 5;
% b = struct('x', 5, 'z', 12, 'y', {b_1 b_1+1});
% c = b;
% %c(1).x = 8;
% d = struct('x', 6, 'y', {b_1});
% 
% eq1 = g*h ==-1;
% eq2 = g+2*h == 3;
% 
% sol = solve(eq1,eq2, [g h]);

%%


% syms x y
% 
% f = (x*y) == -1;
% g = (x)+y == 0;
% solution = testfunct(f, g, x, y)
% 
% 
% function sols1 = testfunct(eq1, eq2, t, s)
% 
%     sols = solve(eq1, eq2, [t s]);
%     var1 = char(t);
%     var2 = char(s);
%     sols1 = sols.(var1)
% %     sols1 = sols.t(1)(imag(sols.t)==0);
% %     sols2 = sols.s(1)(imag(sols.s)==0);
% end
%%
phi = 0:.01:.21;
Positions = phi;
alu = struct('E', 69e9, 'o_adm', 110e6/2);

Materials = [alu];
b = 30e-3;
 pl2 = struct(   'location', 'lev2-lev1',   ...
                 'type',     'lame',        ...
                 'k',        1,             ...
                 'cor_adm',  5,             ...
                 'ener_var', phi,           ... 
                 'Dims',     [3; zeros(7,1)]);
 pl4 = struct(   'location', 'lev1-gnd',    ...
                 'type',     'col',         ...
                 'k',        0.005,             ...
                 'cor_adm',  5,             ...
                 'ener_var', phi,           ... 
                 'Dims',     [b; zeros(7,1)]    );


Pivots = [pl2 pl4];


for i = 1:2;
for j = 1;
switch Pivots(i).type
    case 'lame'
        syms h L
        rig = Pivots(i).k == Pivots(i).Dims(1) * Materials(j).E * b * h^3 / L^3;
        adm = max(abs(Pivots(i).ener_var)) == Materials(j).o_adm*L^2 /(3*Materials(j).E*h);
        fprintf('lame for pivot %d\n',i)
        [h L] = pivotSolve(rig, adm, h, L)
        Pivots(i).dims(1) = h
    case 'col'
        syms e r
        rig = Pivots(i).k == 2* Materials(j).E * b * e^(2.5) / (9*pi*r^(0.5));
        adm = max(abs(Pivots(i).ener_var)) == 3*pi*Materials(j).o_adm*sqrt(r)/(4*Materials(j).E*sqrt(e));
        fprintf('col for pivot %d\n',i)
        [e r] = pivotSolve(rig, adm, e, r)
    end

end
end

function [sol1 sol2] = pivotSolve(eq1, eq2, t, s)
        
        assume(t, {'real', 'positive'});
        assume(s, {'real', 'positive'})
        solutions = vpasolve(eq1, eq2, [t s],'random',true);
        var1 = char(t);
        var2 = char(s);
        
        sol1  = eval(solutions.(var1));
        sol2  = eval(solutions.(var2));
%         sol1_real = eval(solutions.(var1)(imag(solutions.(var1))==0));
%         sol2_real = eval(solutions.(var2)(imag(solutions.(var2))==0));

        %pos_solutions = struct(var1, sol1, var2, sol2);
end

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