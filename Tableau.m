classdef Tableau < handle
properties
    tab
    stab_size
    cir
    wid
    Ns
    boundary
    frozen_qubits
end

methods
    function obj = Tableau(varargin)
        ip = inputParser;
        ip.KeepUnmatched = true;
        ip.addParameter('cir', 2);
        ip.addParameter('wid', 2);
        ip.addParameter('boundary', 'open');
        ip.addParameter('tab', []);
        ip.addParameter('stab_size', 0);
        ip.addParameter('frozen_qubits', []);

        % Parse input arguments
        ip.parse(varargin{:});
        pars = ip.Results;

        % Loop through the fields of pars and assign to properties
        fields = fieldnames(pars);
        for i = 1:numel(fields)
            field = fields{i};
            obj.(field) = pars.(field); % Dynamic field assignment
        end

        obj.Ns = obj.cir * obj.wid * 2;
        if isempty(obj.frozen_qubits)
            obj.frozen_qubits = containers.Map('KeyType', 'int64', 'ValueType', 'any');
        end

        if isempty(obj.tab)
            Ns = obj.Ns;
            tab = eye(Ns*2, Ns*2);
            if strcmp(obj.boundary, 'open')
                tab([1, Ns+1, Ns-1, 2*Ns-1], :) = tab([Ns-1, 2*Ns-1, 1, Ns+1], :);
            end
            obj.tab = tab;
        end
    end
    

    function render_table(obj)
        obj.Ns = obj.cir * obj.wid * 2;
        for i = 1 : obj.Ns*2
            if i == obj.stab_size + 1 
                fprintf('r.....\n')
            end
            if i == obj.Ns + 1
                fprintf('--------\n')
            end
            if i == obj.Ns + obj.stab_size + 1
                fprintf('r.....\n')
            end
            fprintf("%d: %s\n", i, Util.row2pauli(obj.tab(i, :), obj.cir, obj.wid));
        end

        if obj.Ns + obj.stab_size + 1 == obj.Ns*2 + 1
            fprintf('r.....\n')
        end

        fprintf('\n');
    end


    function pair_tab_property(obj)
        obj.Ns = obj.cir * obj.wid * 2;
        b = true;
        str_arr = {};
        for i = 1:obj.Ns*2
            for j = i+1:obj.Ns*2
                if j - obj.Ns == i
                    res = 1 == Util.symplectic_inner_product(obj.tab(i, :), obj.tab(j, :), obj.Ns);
                else
                    res = 0 == Util.symplectic_inner_product(obj.tab(i, :), obj.tab(j, :), obj.Ns);
                end
                if ~res
                    str_arr = [str_arr, [i, j]];
                end
                b = b && res;
            end
        end
        if ~b
            obj.render_table();
            disp(str_arr);
            error('error: not pair satisfy tab property');
        end
    end

end
end

