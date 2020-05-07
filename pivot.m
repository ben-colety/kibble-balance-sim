classdef pivot
    properties
        type
        k
        ener_var
        dims
    end
    methods
        function obj = pivot(type, k, ener_var)
            obj.type = type;
            obj.k = k;
            obj.ener_var = ener_var;
            obj.dims = {};
        end
    end
end