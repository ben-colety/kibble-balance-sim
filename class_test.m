clc
clear all

phi = 0:.01:.1;
Positions = phi;
alu = struct('E', 69e9, 'o_adm', 110e6/2);

Materials = [alu];
b = 50e-3;

pl2 = pivot('spring', 1, phi,1);
pl4 = pivot('point', 1, phi);
pl3 = pl2;

Pivots = [pl2 pl3 pl4];


for i = 1:length(Pivots)
    for j = 1:length(Materials)
        switch Pivots(i).type
            case 'spring'
                syms h L
                rig = Pivots(i).k == Pivots(i).num_lames * Materials(j).E * b * h^3 / L^3;
                adm = max(abs(Pivots(i).ener_var)) == Materials(j).o_adm*L^2 /(3*Materials(j).E*h);
                fprintf('lame for pivot %d\n',i)
                Pivots(i) = pivotSolve(Pivots(i),rig, adm, h, L)

            case 'point'
                %calculation as a col
                syms e r
                rig = Pivots(i).k == 2* Materials(j).E * b * e^(2.5) / (9*pi*r^(0.5));
                adm = max(abs(Pivots(i).ener_var)) == 3*pi*Materials(j).o_adm*sqrt(r)/(4*Materials(j).E*sqrt(e));
                fprintf('col for pivot %d\n',i)
                Pivots(i) = pivotSolve(Pivots(i),rig, adm, e, r)

                %calculation as a cross
                syms h L
                rig2 = Pivots(i).k == 8*Materials(j).E*b*h^3 /(12*L);
                adm2 = max(abs(Pivots(i).ener_var)) == Materials(j).o_adm * L /(2*Materials(j).E*h);
                fprintf('cross for pivot %d\n',i)
                Pivots(i) = pivotSolve(Pivots(i), rig2, adm2, h, L)
        end
    end
end

% function [sol1 sol2] = pivotSolve(eq1, eq2, t, s)
%         assume(t, {'real', 'positive'});
%         assume(s, {'real', 'positive'})
%         solutions = vpasolve(eq1, eq2, [t s],'random',true);
%         var1 = char(t);
%         var2 = char(s);
%         
%         sol1  = eval(solutions.(var1));
%         sol2  = eval(solutions.(var2));
% %         sol1_real = eval(solutions.(var1)(imag(solutions.(var1))==0));
% %         sol2_real = eval(solutions.(var2)(imag(solutions.(var2))==0));
% 
%         %pos_solutions = struct(var1, sol1, var2, sol2);
% end