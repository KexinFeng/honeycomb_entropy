classdef Tracer < handle
properties
    tableau
    cir
    wid
    % Internal
    Ns
end
    
methods
    function obj = Tracer(varargin)
        ip = inputParser;
        ip.addParameter('tableau', []);
        ip.addParameter('cir', 2);
        ip.addParameter('wid', 2);
        
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
        if isempty(obj.tableau)
            obj.tableau = Tableau();
        end
    end
    
    
    function partial_trace(obj, qubits)
        valid_stab_idx = obj.Ns + (1: obj.tableau.stab_size);
        for qubit = qubits
            row_idx_x = valid_stab_idx(obj.tableau.tab(valid_stab_idx, qubit) == 1);
            if ~isempty(row_idx_x)
                obj.add_onto(row_idx_x(1), row_idx_x(2:end));
                valid_stab_idx(valid_stab_idx == row_idx_x(1)) = [];
            end
            row_idx_z = valid_stab_idx(obj.tableau.tab(valid_stab_idx, qubit + obj.Ns) == 1);
            if ~isempty(row_idx_z)
                obj.add_onto(row_idx_z(1), row_idx_z(2:end));
                valid_stab_idx(valid_stab_idx == row_idx_z(1)) = [];
            end
            
            obj.tableau.stab_size = obj.tableau.stab_size - ~isempty(row_idx_x) - ~isempty(row_idx_z);
        end
        assert(obj.tableau.stab_size == length(valid_stab_idx));
        
        src = obj.Ns + (1: obj.tableau.stab_size);
        obj.tableau.tab([src, valid_stab_idx], :) = ...
            obj.tableau.tab([valid_stab_idx, src], :);
        obj.tableau.tab([src - obj.Ns, valid_stab_idx - obj.Ns], :) = ...
            obj.tableau.tab([valid_stab_idx - obj.Ns, src - obj.Ns], :);
    end


    function add_onto(obj, row_src, row_tgt)
        if isempty(row_tgt)
            return
        end
        row = obj.tableau.tab(row_src, :);
        obj.tableau.tab(row_tgt, :) = mod(obj.tableau.tab(row_tgt, :) + row, 2);
        
        % destab update
        row_src_bar = row_src - obj.Ns;
        row_tgt_bar = row_tgt - obj.Ns;
        obj.tableau.tab(row_src_bar, :) = mod(obj.tableau.tab(row_src_bar, :) + sum(obj.tableau.tab(row_tgt_bar, :), 1), 2);
    end
end
end

