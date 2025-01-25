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
        ip.PartialMatching = false;
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
            tab = uint16(eye(Ns*2, Ns*2));
            if canUseGPU()
                tab = gpuArray(tab);
            end
            if strcmp(obj.boundary, 'open')
            % Distinction relies on stab_size, if =N_s pure, if =0 mixed.
                tab([1, Ns+1, Ns-1, 2*Ns-1], :) = tab([Ns-1, 2*Ns-1, 1, Ns+1], :);
                obj.frozen_qubits(1) = true;
                obj.frozen_qubits(Ns+1) = true;
            end
            obj.tab = tab;
        end
    end


    function tableau = clone(obj)
        kv_list = get_param(obj);

        if isempty(obj.frozen_qubits)
            clonedMap = containers.Map(); 
        else
            clonedMap = containers.Map(keys(obj.frozen_qubits), values(obj.frozen_qubits));
        end
        % This should reproduce the same behavior, yes? Tested with both
        % empty and occupied frozen_qubits map and seemed to work fine
        tableau = Tableau(kv_list{:});%, 'frozen_qubits', clonedMap);
    end
    

    function save(file_path)
        kv_list = get_param(obj);
        save(file_path, 'kv_list');
    end

    function entropy = get_entropy(obj)
        % entropy = obj.Ns - 2*strcmp(obj.boundary, 'open') - obj.stab_size; % frozen qubits are discarded.
        entropy = obj.Ns - length(keys(obj.frozen_qubits)) - obj.stab_size; % frozen qubits are discarded.
    end
    
    %% Tracer
    function partial_trace(obj, qubits)
        % qubits: [1, |qubits|]
        
        % Performs trace over sites with indices given in [qubits]
        Ns = obj.Ns;
        stab_size = obj.stab_size;

        % search
        valid_stab_idx = Ns + (1: stab_size); %array of all stab row indices
        disp(valid_stab_idx);
        for qubit = qubits
            if isKey(obj.frozen_qubits, qubit)
                % Does not trace out frozen qubits
                continue
            end
            % Freezes traced qubits
            obj.frozen_qubits(qubit) = true;

            % Find all Xs or Ys on the qubit, transform them and remove the first
            row_idx_x = valid_stab_idx(obj.tab(valid_stab_idx, qubit) == 1);
            disp({'x';row_idx_x});
            if ~isempty(row_idx_x)
                % If gate from current qubit is present in another row,
                % removes that influence
                obj.add_onto(row_idx_x(1), row_idx_x(2:end));
                % Removes the traced qubit from list
                disp(row_idx_x(1));
                valid_stab_idx(valid_stab_idx == row_idx_x(1)) = [];
                disp(size(valid_stab_idx));
                disp(valid_stab_idx);
            end
            % Find all Zs on the qubit, transform them and remove the first
            row_idx_z = valid_stab_idx(obj.tab(valid_stab_idx, qubit + Ns) == 1);
            disp({'z';row_idx_z});
            if ~isempty(row_idx_z)
                obj.add_onto(row_idx_z(1), row_idx_z(2:end));
                disp(row_idx_z(1));
                valid_stab_idx(valid_stab_idx == row_idx_z(1)) = [];
                disp(size(valid_stab_idx));
                disp(valid_stab_idx);
            end
            disp({row_idx_x;row_idx_z});
            disp([~isempty(row_idx_x) ~isempty(row_idx_z)]);
            % Due to y gates being sorted by row_idx_x, this only ever
            % subtracts 1 from the stab size, as we expect.
            stab_size = stab_size - ~isempty(row_idx_x) - ~isempty(row_idx_z);
            disp(stab_size);
        end
        assert(stab_size == length(valid_stab_idx));
        
        % Re-org to remove the X and Z row from the stabilizer section
        % Reorder tableau with only the untraced stabilizers beginning at
        % row N_s+1. Traced columns move to the end in both stab and
        % destabilizer sections, essentially, not totally clear on
        % intended effect here
        src = Ns + (1: obj.stab_size);
        new_order = [valid_stab_idx, setdiff(src, valid_stab_idx)];
        disp(new_order);
        % Is this not how the rows should be reassigned? Testing it this
        % way seemed to put the expected qubits outside the stabilizer.
        %obj.tab(new_order, :) = obj.tab(src, :);        
        %obj.tab(new_order - Ns, :) = obj.tab(src - Ns, :);
        obj.tab(src, :) = obj.tab(new_order, :) ;        
        obj.tab(src - Ns, :) = obj.tab(new_order - Ns, :);
        obj.stab_size = stab_size;
    end


    function add_onto(obj, row_src, row_tgt)
        % Adds row_src onto row_tgt and then adds the corresponding
        % destabilizers together onto the row_src complement
        % works for multiple row_tgt simulatenously
        if isempty(row_tgt)
            return
        end
        Ns = obj.Ns;

        row = obj.tab(row_src, :);
        obj.tab(row_tgt, :) = mod(obj.tab(row_tgt, :) + row, 2);
        
        % destab update
        row_src_bar = row_src - Ns;
        row_tgt_bar = row_tgt - Ns;
        obj.tab(row_src_bar, :) = mod(obj.tab(row_src_bar, :) + uint16(sum(obj.tab(row_tgt_bar, :), 1)), 2);
    end
   
    %% Utility
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
        %This is a duplicate of the Util.m function? Checks that one of the
        %model assumptions is satisfied: only rows i and i+Ns anticommute.
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

