classdef GaussianEliminator < handle
properties
    tab
    stab_size
    cir
    wid
    % Internal
    Ns
end
    
methods
    function obj = GaussianEliminator(varargin)
        ip = inputParser;
        ip.addParameter('tab', []);
        ip.addParameter('stab_size', []);
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
    end
    
    
    function partial_trace(obj, qubits)
        valid_stab_idx = obj.Ns + (1: obj.stab_size);
        for qubit = qubits
            row_idx_x = valid_stab_idx(obj.tab(valid_stab_idx, qubit) == 1);
            if ~isempty(row_idx_x)
                obj.add_onto(row_idx_x(1), row_idx_x(2:end));
                valid_stab_idx(valid_stab_idx == row_idx_x(1)) = [];
            end
            row_idx_z = valid_stab_idx(obj.tab(valid_stab_idx, qubit + obj.Ns) == 1);
            if ~isempty(row_idx_z)
                obj.add_onto(row_idx_z(1), row_idx_z(2:end));
                valid_stab_idx(valid_stab_idx == row_idx_z(1)) = [];
            end
            
            obj.stab_size = obj.stab_size - ~isempty(row_idx_x) - ~isempty(row_idx_z);
        end
        assert(obj.stab_size == length(valid_stab_idx));
        
        src = obj.Ns + (1: obj.stab_size);
        obj.tab([src, valid_stab_idx], :) = ...
            obj.tab([valid_stab_idx, src], :);
        obj.tab([src - obj.Ns, valid_stab_idx - obj.Ns], :) = ...
            obj.tab([valid_stab_idx - obj.Ns, src - obj.Ns], :);
    end


    function add_onto(obj, row_src, row_tgt)
        if isempty(row_tgt)
            return
        end
        row = obj.tab(row_src, :);
        obj.tab(row_tgt, :) = mod(obj.tab(row_tgt, :) + row, 2);
        
        % destab update
        row_src_bar = row_src - obj.Ns;
        row_tgt_bar = row_tgt - obj.Ns;
        obj.tab(row_src_bar, :) = mod(obj.tab(row_src_bar, :) + sum(obj.tab(row_tgt_bar, :), 1), 2);
    end
end
end

