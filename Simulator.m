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
        
        obj.tableau = Tableau(varargin{:});

        if obj.verbose
            obj.tableau.render_table();
        end

        obj.check_generator = CheckGenerator(varargin{:});
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

        % Get the list of properties
        propNames = properties(obj);
        numProps = length(propNames);
        kvList = cell(1, 2 * numProps); 
        
        % Construct the key-value list
        for i = 1:numProps
            kvList{2*i - 1} = propNames{i};  
            kvList{2*i} = obj.(propNames{i});
        end

        obj.check_generator = CheckGenerator(kvList{:});
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
        
        tableau = obj.tableau;
        % tab = obj.tableau.tab;
        % stab_size = obj.tableau.stab_size;

        %% main loop
        for step = 1:T
            for idx = 1:cir*wid %#ok<*PROP>
                [y, x] = homod(idx, wid);
                % [tableau, stab_size, bond] = measure_rank(x, y, tableau, stab_size, pars);
                bond = obj.measure_destab(x, y);
                if isempty(bond)
                    continue
                end
                
                if obj.verbose
                    % Util.pair_tab_property(tab, cir, wid);
                    tableau.pair_tab_property();
                    % print
                    fprintf('\nidx: %d\n', idx);
                    tableau.render_table();
                end
            end
    
            if obj.verbose
                % Util.pair_tab_property(tab, cir, wid);
                tableau.pair_tab_property()
                fprintf('step=%d, r=%d\n\n', step, tableau.stab_size)
                disp(' ')
            end
    
            % Plotting
            if obj.plotting
                plot(step, Ns - 2*strcmp(obj.boundary, 'open') - obj.tableau.stab_size, '*', 'Color', 'b');
                pause(0.01);
            end
        end
        
        disp(['stab_size=', num2str(tableau.stab_size), ' total spin:', num2str(Ns - 2*strcmp(obj.boundary, 'open'))])
        
        % Save
        tab = tableau.tab;
        stab_size = tableau.stab_size;
        folder = sprintf('./tmp/');
        mkdir(folder);
        save(strjoin({folder, 'tab4unit_test.mat'}, ""), "tab", "stab_size")
    end


    function bond = measure_destab(obj, x, y)
        % Measure on the unit cell (x, y)     
        cir = obj.cir;
        wid = obj.wid;
    
        
        [row, bond] = obj.check_generator.generate_bond(x, y);
        if isempty(bond)
            return
        end
    
        if obj.verbose
           fprintf('bond: %s\n', Util.row2pauli(row, cir, wid))
        end
    
        
        [scenario, row_idx] = obj.check_scenario(row);
        
        if scenario == 1
            % stochastic output containing de-stablizer
            obj.scenario1(row, row_idx);
        elseif scenario == 2
            % deterministic output
            % ;
        else
            % stochastic output containing enhanced space
            obj.scenario3(row, row_idx);
        end
     
    end
    
    
    %% Tableaue processing
    function scenario1(obj, row_measure, row_idx)
        tableau = obj.tableau;        
        % row_idx points to the row anticommuting with row_measure
        Ns = size(row_measure, 2) / 2;
        
        row_append = tableau.tab(row_idx, :);
        % restore the tableau property
        for i = [row_idx + 1: Ns*2, 1: Ns]
            if Util.symplectic_inner_product(tableau.tab(i, :), row_measure, Ns)
                tableau.tab(i, :) = Util.pauli_product(row_append, tableau.tab(i, :));
            end
        end
        tableau.tab(row_idx, :) = row_measure;
        tableau.tab(row_idx - Ns, :) = row_append;
    end
    

    function scenario3(obj, row_measure, row_idx)
        tableau = obj.tableau;
        stab_size = tableau.stab_size;
        
        % row_idx points to the row anticommuting with row_measure
        Ns = size(row_measure, 2) / 2;
    
        % swap row_idx to Ns + stab_size + 1
        row_idx_bar = homod(row_idx + Ns, 2*Ns);
        tableau.tab([row_idx, Ns + stab_size + 1], :) = tableau.tab([Ns + stab_size + 1, row_idx], :);
        if row_idx_bar ~= Ns + stab_size + 1
            tableau.tab([row_idx_bar, stab_size + 1], :) = tableau.tab([stab_size + 1, row_idx_bar], :);
        end
        row_idx = Ns + stab_size + 1;
    
        obj.scenario1(row_measure, row_idx);
        tableau.stab_size = stab_size + 1;
    end
    
    
    function [scenario, row_idx] = check_scenario(obj, row_measure)
        Ns = obj.tableau.Ns;
        tableau = obj.tableau;
        stab_size = obj.tableau.stab_size;
        
        Nrow = Ns - triexp(strcmp(obj.boundary, 'open'), 2, 0);
    
        % scenario 1
        for i = Ns+1: Ns + stab_size
            if Util.symplectic_inner_product(tableau.tab(i, :), row_measure, Ns)
                scenario = 1;
                row_idx = i;
                return
            end
        end
        
        % scenario 3
        for i = [Ns + stab_size + 1: Ns + Nrow, stab_size + 1: Nrow]
            if Util.symplectic_inner_product(tableau.tab(i, :), row_measure, Ns)
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
    
end
end

