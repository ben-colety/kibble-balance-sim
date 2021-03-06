format shorteng
phi = 0:.01:.12;
Positions = phi;
alu = struct('E', 69e9, 'o_adm', 110e6/2);

Materials = [alu];
b = 50e-3;

pl2 = pivot('spring', 1, phi,1,1);
pl4 = pivot('point', 4, phi, 1);
pc1 = pivot('parallel', 5, phi, 1);
pl3 = pl2;

Pivots = [pl2 pl3 pl4 pc1];


for i = 1:length(Pivots)
    for j = 1:length(Materials)
        switch Pivots(i).type
            case {'spring','parallel'}
                syms h L
                cond1 = L <= 0.1;
                rig = Pivots(i).k == Pivots(i).num_lames * Materials(j).E * b * h^3 / L^3;
                adm = max(abs(Pivots(i).ener_var)) == Materials(j).o_adm*L^2 /(3*Materials(j).E*h);
                if Pivots(i).type == "parallel"
                    fprintf('parallel for pivot %d\n',i)
                elseif Pivots(i).type == "spring"
                    fprintf('spring for pivot %d with %d lames\n',i,Pivots(i).num_lames)
                end
                Pivots(i) = pivotSolve(Pivots(i),rig, adm, h, L,cond1)

            case {'point','col','cross'}
                %calculation as a col
                syms e r
                cond2 = r <= 0.01;
                rig = Pivots(i).k == 2* Materials(j).E * b * e^(2.5) / (9*pi*r^(0.5));
                adm = max(abs(Pivots(i).ener_var)) == 3*pi*Materials(j).o_adm*sqrt(r)/(4*Materials(j).E*sqrt(e));
                fprintf('col for pivot %d\n',i)
                Pivots(i) = pivotSolve(Pivots(i),rig, adm, e, r,cond2)

                %calculation as a cross
                syms h L
                cond3 = L <= 60*h;
                rig2 = Pivots(i).k == 8*Materials(j).E*b*h^3 /(12*L);
                adm2 = max(abs(Pivots(i).ener_var)) == Materials(j).o_adm * L /(2*Materials(j).E*h);
                fprintf('cross for pivot %d\n',i)
                Pivots(i) = pivotSolve(Pivots(i), rig2, adm2, h, L,cond3)
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