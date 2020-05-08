classdef pivot
    properties
        type;
        k;
        ener_var;
        num_lames = 0;
        quant = 1;
        h; L; e; r;
        
    end
    methods
        function obj = pivot(type, k, ener_var, quant, num_lames)
            obj.type = type;
            obj.k = k;
            obj.ener_var = ener_var;
            if quant > 1
                obj.quant = quant;
            end
            if type == "spring"
                obj.num_lames = num_lames;
            elseif type == "parallel"
                obj.num_lames = 2;
            end
        end
        function obj = pivotSolve(obj, eq1, eq2, t, s)
            assume(t, {'real', 'positive'});
            assume(s, {'real', 'positive'})
            solutions = solve(eq1, eq2, [t s]);
            var1 = char(t);
            var2 = char(s);

            sol1  = eval(solutions.(var1));
            sol2  = eval(solutions.(var2));
            obj.(var1) = sol1;
            obj.(var2) = sol2;
            %obj.dims = {sol1 sol2}
        end
        function energy = calcEnergy(p)
            if class(p.ener_var) == "sym"
                energy = 0.5*p.k*p.ener_var^2;
            elseif class(p.ener_var) == "double"
                energy = 0.5*p.k*p.ener_var.^2;
            end
        end
    end
end