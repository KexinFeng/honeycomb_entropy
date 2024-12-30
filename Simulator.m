classdef Simulator < handle
properties
    cir
    wid
    T
    plotting
    boundary
    shift
    verbose
    %% internal 
    tableau
    check_generator
end

methods
    function obj = Simulator(varargin)
        ip = inputParser;
        ip.KeepUnmatched = true;
        ip.PartialMatching = false;
        ip.addParameter('cir', 10);
        ip.addParameter('wid', 10);
        ip.addParameter('boundary', 'open');
        ip.addParameter('shift', 0);
        ip.addParameter('T', 30);
        ip.addParameter('plotting', 1);
        ip.addParameter('verbose', false);
    
        ip.parse(varargin{:});
        pars = ip.Results;
        
        pars.wid = pars.cir;

        % Loop through the fields of pars and assign to properties
        fields = fieldnames(pars);
        for i = 1:numel(fields)
            field = fields{i};
            obj.(field) = pars.(field);  % Dynamic field assignment
        end

        % Convert pars to a cell array of parameter-value pairs
        values = struct2cell(pars);
        fields2 = fieldnames(ip.Unmatched);
        values2 = struct2cell(ip.Unmatched);
        varargin = reshape([[fields; fields2], [values; values2]]', 1, []); 

        % Internal
        obj.tableau = Tableau(varargin{:});
        if obj.verbose
            obj.tableau.render_table();
        end

        obj.check_generator = CheckGeneratorPlaq(varargin{:});
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

        kvList = get_param(obj);
        obj.check_generator = CheckGenerator(kvList{:});
    end
    

    function [ys, xs] = simulate(obj)
        tic 
        filepath = fileparts(mfilename('fullpath'));
        addpath(genpath([filepath, '/../utils']));
        % addpath([filepath, '/util.m']);

    
        %% prepare
        cir = obj.cir;
        wid = obj.wid;
        Ns = cir*wid*2;
        
        T = obj.T;

        xs = zeros(1, T);
        ys = zeros(1, T);
        
        if obj.plotting
            prep_plots;
            figure
            hold on

            ylim([0, inf]);
            xlim([0, obj.T]);

            plot(0, Ns - 2*strcmp(obj.boundary, 'open'), '*', 'Color', 'b');
            pause(0.01);        
        end
        
        tableau = obj.tableau;

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
                    fprintf('\nidx: %d\n', idx);
                    % tableau.render_table();
                end
            end
    
            if obj.verbose
                % Util.pair_tab_property(tab, cir, wid);
                tableau.pair_tab_property()
                fprintf('step=%d, stab_size=%d\n\n', step, tableau.stab_size)
                disp(' ')
            end
    
            % Plotting
            if obj.plotting
                plot(step, obj.tableau.get_entropy(), '*', 'Color', 'b');
                pause(0.01);
            end

            % Output
            ys(step) = obj.tableau.get_entropy();
            xs(step) = step;
        end
        
        tableau.render_table()
        disp(['stab_size=', num2str(tableau.stab_size), ' total spin:', num2str(Ns - 2*strcmp(obj.boundary, 'open'))])
        disp(['entropy=', num2str(tableau.get_entropy())])
        
        % % Save
        % tab = tableau.tab;
        % stab_size = tableau.stab_size;
        % folder = sprintf('./tmp/');
        % mkdir(folder);
        % save(strjoin({folder, sprintf('tmp%d.mat', obj.T)}, ""), "tab", "stab_size")
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
           fprintf('check: %s\n', Util.row2pauli(row, cir, wid))
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
        
        % Find the first non-commuting row
        row_append = tableau.tab(row_idx, :);
        % Restore the tableau property for the rest non-commuting rows
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
    
        % Swap row_idx to Ns + stab_size + 1
        tableau.tab([row_idx, Ns + stab_size + 1], :) = tableau.tab([Ns + stab_size + 1, row_idx], :);
        % Swap row_idx_bar accordingly
        row_idx_bar = homod(row_idx + Ns, 2*Ns);
        if row_idx_bar ~= Ns + stab_size + 1
            % Deduplicate the swap
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

end
end

