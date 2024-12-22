classdef Simulator < handle
properties
    cir
    wid
    plotting
    boundary
    verbose
    T
    %% internal 
    tableau
    check_generator
end

methods
    function obj = Simulator(varargin)
        ip = inputParser;
        ip.addParameter('cir', 10);
        ip.addParameter('wid', 10);
        ip.addParameter('boundary', 'open');
        ip.addParameter('T', 30);
        ip.addParameter('plotting', 1);
        ip.addParameter('verbose', false);
    
        ip.parse(varargin{:});
        pars = ip.Results;

        % Loop through the fields of pars and assign to properties
        fields = fieldnames(pars);
        for i = 1:numel(fields)
            field = fields{i};
            obj.(field) = pars.(field);  % Dynamic field assignment
        end
        
        obj.tableau = Tableau('cir', obj.cir, 'wid', obj.wid, 'boundary', obj.boundary);

        if obj.verbose
            obj.tableau.render_table();
        end

        obj.check_generator = [];
    end


    function set_tableau(obj, tableau)
        % Loop through the fields of pars and assign to properties
        fields = fieldnames(tableau);
        for i = 1:numel(fields)
            field = fields{i};
            if isprop(obj, field) % Check if the field exists as a property in the object
                obj.(field) = tableau.(field);  % Dynamic field assignment
            end
        end
        obj.tableau = tableau;
    end
    

    function simulate(obj)
        tic 
        filepath = fileparts(mfilename('fullpath'));
        addpath(genpath([filepath, '/../utils']));
        % addpath([filepath, '/util.m']);
        
        rng(24);
        clc;
    
        %% prepare
        cir = obj.cir;
        wid = obj.wid;
        Ns = cir*wid*2;
        
        T = obj.T;
        
        if obj.plotting
            figure
            hold on
            plot(0, Ns - 2*strcmp(obj.boundary, 'open'), '*', 'Color', 'b');
            pause(0.01);        
        end
        
        tableau = obj.tableau.tab;
        stab_size = obj.tableau.stab_size;

        %% main loop
        for step = 1:T
            for idx = 1:cir*wid %#ok<*PROP>
                [y, x] = homod(idx, wid);
                % [tableau, stab_size, bond] = measure_rank(x, y, tableau, stab_size, pars);
                [tableau, stab_size, bond] = obj.measure_destab(x, y, tableau, stab_size);
                if isempty(bond)
                    continue
                end
                
                if obj.verbose
                    Util.pair_tab_property(tableau, cir, wid);
                    % print
                    fprintf('\nidx: %d\n', idx);
                    Util.render_table_destab(tableau, stab_size, cir, wid);
                end
            end
    
            if obj.verbose
                Util.pair_tab_property(tableau, cir, wid);
                fprintf('step=%d, r=%d\n\n', step, stab_size)
                disp(' ')
            end
    
            % Plotting
            if obj.plotting
                plot(step, Ns - 2*strcmp(obj.boundary, 'open') - stab_size, '*', 'Color', 'b');
                pause(0.01);
            end
        end
        
        disp(['stab_size=', num2str(stab_size), ' total spin:', num2str(Ns - 2*strcmp(obj.boundary, 'open'))])
        
        % Save
        folder = sprintf('./tmp/');
        mkdir(folder);
        save(strjoin({folder, 'tab4unit_test.mat'}, ""), "tableau", "stab_size")
    end


    function [tab, stab_size, bond] = measure_destab(obj, x, y, tab, stab_size)
        % Measure on the unit cell (x, y)     
        cir = obj.cir;
        wid = obj.wid;
        Ns = cir*wid*2;
    
    
        [row, bond] = obj.generate_bond(x, y);
        if isempty(bond)
            return
        end
    
        if obj.verbose
           fprintf('bond: %s\n', Util.row2pauli(row, cir, wid))
        end
    
        
        [scenario, row_idx] = obj.check_scenario(tab, row, stab_size);
        
        if scenario == 1
            % stochastic output containing de-stablizer
            tab = obj.scenario1(tab, row, row_idx);
        elseif scenario == 2
            % deterministic output
            % ;
        else
            % stochastic output containing enhanced space
            [tab, stab_size] = obj.scenario3(tab, row, row_idx, stab_size);
        end
     
    end
    
    
    %% Tableaue processing
    function tab = scenario1(~, tab, row_measure, row_idx)
        % row_idx points to the row anticommuting with row_measure
        Ns = size(row_measure, 2) / 2;
        
        row_append = tab(row_idx, :);
        % restore the tableau property
        for i = [row_idx + 1: Ns*2, 1: Ns]
            if Util.symplectic_inner_product(tab(i, :), row_measure, Ns)
                tab(i, :) = Util.pauli_product(row_append, tab(i, :));
            end
        end
        tab(row_idx, :) = row_measure;
        tab(row_idx - Ns, :) = row_append;
    end
    

    function [tab, stab_size] = scenario3(obj, tab, row_measure, row_idx, stab_size)
        % row_idx points to the row anticommuting with row_measure
        Ns = size(row_measure, 2) / 2;
    
        % swap row_idx to Ns + stab_size + 1
        row_idx_bar = homod(row_idx + Ns, 2*Ns);
        tab([row_idx, Ns + stab_size + 1], :) = tab([Ns + stab_size + 1, row_idx], :);
        if row_idx_bar ~= Ns + stab_size + 1
            tab([row_idx_bar, stab_size + 1], :) = tab([stab_size + 1, row_idx_bar], :);
        end
        row_idx = Ns + stab_size + 1;
    
        tab = obj.scenario1(tab, row_measure, row_idx);
        stab_size = stab_size + 1;
    end
    
    
    function [scenario, row_idx] = check_scenario(obj, tab, row_measure, stab_size)
        Ns = obj.cir * obj.wid * 2;
        
        Nrow = Ns - triexp(strcmp(obj.boundary, 'open'), 2, 0);
    
        % scenario 1
        for i = Ns+1: Ns + stab_size
            if Util.symplectic_inner_product(tab(i, :), row_measure, Ns)
                scenario = 1;
                row_idx = i;
                return
            end
        end
        
        % scenario 3
        for i = [Ns + stab_size + 1: Ns + Nrow, stab_size + 1: Nrow]
            if Util.symplectic_inner_product(tab(i, :), row_measure, Ns)
                scenario = 3;
                row_idx = i;
                return
            end
        end
    
        % scenario 2
        scenario = 2;
        row_idx = 0;
    end

    %% Tracer
    function partial_trace(obj, qubits)
        tableau = obj.tableau;
        Ns = obj.tableau.Ns;
        stab_size = obj.tableau.stab_size;

        % search
        valid_stab_idx = Ns + (1: stab_size);
        for qubit = qubits
            row_idx_x = valid_stab_idx(tableau.tab(valid_stab_idx, qubit) == 1);
            if ~isempty(row_idx_x)
                obj.add_onto(row_idx_x(1), row_idx_x(2:end));
                valid_stab_idx(valid_stab_idx == row_idx_x(1)) = [];
            end
            row_idx_z = valid_stab_idx(tableau.tab(valid_stab_idx, qubit + Ns) == 1);
            if ~isempty(row_idx_z)
                obj.add_onto(row_idx_z(1), row_idx_z(2:end));
                valid_stab_idx(valid_stab_idx == row_idx_z(1)) = [];
            end
            
            stab_size = stab_size - ~isempty(row_idx_x) - ~isempty(row_idx_z);
        end
        assert(stab_size == length(valid_stab_idx));
        
        src = Ns + (1: stab_size);
        tableau.tab([src, valid_stab_idx], :) = ...
            tableau.tab([valid_stab_idx, src], :);
        tableau.tab([src - Ns, valid_stab_idx - Ns], :) = ...
            tableau.tab([valid_stab_idx - Ns, src - Ns], :);

        % obj.tableau = tableau;
        tableau.stab_size = stab_size;
    end


    function add_onto(obj, row_src, row_tgt)
        if isempty(row_tgt)
            return
        end
        tableau = obj.tableau;
        Ns = obj.tableau.Ns;

        row = tableau.tab(row_src, :);
        tableau.tab(row_tgt, :) = mod(tableau.tab(row_tgt, :) + row, 2);
        
        % destab update
        row_src_bar = row_src - Ns;
        row_tgt_bar = row_tgt - Ns;
        tableau.tab(row_src_bar, :) = mod(tableau.tab(row_src_bar, :) + sum(tableau.tab(row_tgt_bar, :), 1), 2);
    
        % obj.tableau = tableau;
    end
    

    %% Random generation
    function [row, bond] = generate_bond(obj, x, y)
        cir = obj.cir;
        wid = obj.wid;
        
        % Generate bond
        if strcmp(obj.boundary, "periodic")
            bond = randi([1, 3]); % 1:z, 2:x, 3:y
        elseif strcmp(obj.boundary, 'open')
            if x==1 && y==1 
                bond = randsample([1, 3], 1);
            elseif x == cir && y == wid
                bond = [];
                row = [];
                return
            elseif y == wid
                bond = randsample([2, 3], 1);
            elseif x == cir
                bond = randsample([1, 2], 1);
            else
                bond = randi([1, 3]);
            end
        else
            error('%s not right', obj.boundary)
        end
        
        row = obj.bond2row(x, y, bond);
    end


    function row = bond2row(obj, x, y, bond)
        % bond_names = ['Z', 'X', 'Y'];
    
        cir = obj.cir;
        wid = obj.wid;
        Ns = cir * wid * 2;
    
        row = zeros(1, Ns*2);
        % Site index: x = 1~cir; y = 1~wid; ab = 0,1
        pos = @(x, y, ab) y + wid.*(x-1) + cir*wid.*ab;
    
        % (x, y, 0)-X-(x, y , 1)
        % (x, y+1, 0)-Z-(x, y, 1)
        % (x+1, y, 0)-Y-(x, y, 1)
        bin = dec2bin(bond, 2);
        bin_arr = kron(str2num(bin(:)), [1; 1]);
        
        xs = [x, homod(x+(bond == 3), cir)];
        ys = [y, homod(y+(bond == 1), wid)];
        ab = [1, 0];
        
        idx = pos(xs, ys, ab);
        idx = [idx, idx + Ns];
        
        row(idx) = reshape(bin_arr, size(idx));
        % bond_names = ['Z', 'X', 'Y'];
        % pauli = Util.row2pauli(row, cir, wid);
        % fprintf('%s\n pauli: %s\n', bond_names{bond}, pauli)
    end

end
end

