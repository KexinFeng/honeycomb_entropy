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
    


    function [es, ts] = zero_flux_init(obj)
        cir = obj.cir;
        wid = obj.wid;

        for idx = 1: cir * wid - 1
            [y, x] = homod(idx, wid);
            row_plaq = obj.check_generator.plaq2row_commutive(x, y);
            [scenario, row_idx] = obj.check_scenario(row_plaq);
            assert(scenario == 3);
            
            obj.scenario3(row_plaq, row_idx);
        end
        
        if obj.verbose
            obj.tableau.render_table();
        end

        ts = 0;
        es = obj.tableau.get_entropy();
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
        % addpath('/Users/kx/Desktop/forked/MatlabProgressBar/');
        % pb = ProgressBar(T, 'IsParallel', false, ...
        %     'Title', sprintf('%d processing...', obj.cir));

        for step = 1: T
            for idx = 1: cir*wid %#ok<*PROP>
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
            
            % Progress
            % pb([], [], []);
            % updateParallel([], pwd);
            fprintf('%d / %d progress %d\n', step, T, obj.cir);
        end

        % pb.release();
        % close(h);

        if obj.verbose
            tableau.render_table()
        end
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

        if canUseGPU()
            row = gpuArray(row);
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
        tab = obj.tableau.tab;
        % tab2 = tab;
        % row_idx points to the row anticommuting with row_measure
        Ns = size(row_measure, 2) / 2;

        % Find the first non-commuting row
        row_append = tab(row_idx, :);

        % % Restore the tableau property for the rest non-commuting rows
        % for i = [row_idx + 1: Ns*2, 1: Ns] 
        %     if Util.symplectic_inner_product(tab(i, :), row_measure, Ns)
        %         tab(i, :) = Util.pauli_product(row_append, tab(i, :));
        %     end
        % end

        % % Restore the tableau property for the rest non-commuting rows
        % % Define the range of indices
        % indices = [row_idx + 1: Ns * 2, 1: Ns];
        % % Extract the relevant rows
        % rows_to_check = tab(indices, :);
        % % Compute the symplectic inner product for all rows in one go
        % phases = Util.symplectic_inner_product_vec(rows_to_check, row_measure, Ns);
        % % Logical mask for rows where symplectic_inner_product equals 1
        % rows_to_update = phases == 1;
        % % Apply the pauli_product operation to the selected rows
        % tab(indices(rows_to_update), :) = mod(row_append + rows_to_check(rows_to_update, :), 2);
      
        % Restore the tableau property for the rest non-commuting rows
        % Define the range of indices
        indices_bool = false(Ns*2, 1);
        indices_bool([row_idx + 1: Ns * 2, 1: Ns]) = true;
        % Compute the symplectic inner product for all rows in one go
        phases = Util.symplectic_inner_product_vec(tab, row_measure, Ns);
        % Logical mask for rows where symplectic_inner_product equals 1
        rows_to_update = (phases == 1 & indices_bool);
        % Apply the pauli_product operation to the selected rows
        tab(rows_to_update, :) = mod(row_append + tab(rows_to_update, :), 2);
        
        % assert(all(tab == tab2, "all"))

        % Assign the updated rows
        tab(row_idx, :) = row_measure;
        tab(row_idx - Ns, :) = row_append;

        obj.tableau.tab = tab;
    end


    function scenario3(obj, row_measure, row_idx)
        tab = obj.tableau.tab;
        stab_size = obj.tableau.stab_size;

        % row_idx points to the row anticommuting with row_measure
        Ns = size(row_measure, 2) / 2;

        % Swap row_idx to Ns + stab_size + 1
        tab([row_idx, Ns + stab_size + 1], :) = tab([Ns + stab_size + 1, row_idx], :);
        % Swap row_idx_bar accordingly
        row_idx_bar = homod(row_idx + Ns, 2*Ns);
        if row_idx_bar ~= Ns + stab_size + 1
            % Deduplicate the swap
            tab([row_idx_bar, stab_size + 1], :) = tab([stab_size + 1, row_idx_bar], :);
        end
        row_idx = Ns + stab_size + 1;
        obj.tableau.tab = tab;

        obj.scenario1(row_measure, row_idx);
        obj.tableau.stab_size = stab_size + 1;
    end

    
    function [scenario, row_idx] = check_scenario(obj, row_measure)
        Ns = obj.tableau.Ns;
        tableau = obj.tableau;
        stab_size = obj.tableau.stab_size;
<<<<<<< Updated upstream

=======
        %What is triexp?
>>>>>>> Stashed changes
        Nrow = Ns - triexp(strcmp(obj.boundary, 'open'), 2, 0);

        phases = Util.symplectic_inner_product_vec(tableau.tab, row_measure, Ns);
        

        % Scenario1
        indices_bool_s1 = false(2*Ns, 1);
        if canUseGPU()
            indices_bool_s1 = gpuArray(indices_bool_s1);
        end
        indices_bool_s1(Ns+1 : Ns+stab_size) = true;
  
        % check if any row satisfies the condition
        row_idx_s1 = find(phases == 1 & indices_bool_s1, 1);
        if ~isempty(row_idx_s1)
            scenario = 1;
            row_idx = row_idx_s1;
            return;
        end
        
        % Scenario3
        indices_bool_s3 = false(2*Ns, 1);
        if canUseGPU()
            indices_bool_s3 = gpuArray(indices_bool_s3);
        end
        indices_bool_s3([Ns + stab_size + 1 : Ns + Nrow, stab_size + 1 : Nrow]) = true;
       
        % check if any row satisfies the condition
        check_bool = phases == 1 & indices_bool_s3;
        % flip the two sections of indices
        check_bool([1: Ns, Ns + 1: 2*Ns]) = check_bool([Ns + 1: 2*Ns, 1: Ns]);
        row_idx_s3 = find(check_bool, 1);
        if ~isempty(row_idx_s3)
            scenario = 3;
            % restore the index by computing its conjugate
            row_idx = mod(row_idx_s3 - 1 + Ns, 2*Ns) + 1;
            return;
        end

        % scenario 2
        scenario = 2;
        row_idx = 0;
    end

    
    function [scenario, row_idx] = check_scenario_0(obj, row_measure)
        Ns = obj.tableau.Ns;
        tableau = obj.tableau;
        stab_size = obj.tableau.stab_size;

        Nrow = Ns - triexp(strcmp(obj.boundary, 'open'), 2, 0);

        % Scenario1
        rows_scenario1 = tableau.tab(Ns+1 : Ns+stab_size, :);

        % Compute the symplectic inner product for all rows in scenario 1
        phases_scenario1 = Util.symplectic_inner_product_vec(rows_scenario1, row_measure, Ns);

        % Check if any row satisfies the condition
        idx_scenario1 = find(phases_scenario1 == 1, 1);
        if ~isempty(idx_scenario1)
            scenario = 1;
            row_idx = Ns + idx_scenario1;
            return;
        end

        % Scenario3
        indices_scenario3 = [Ns + stab_size + 1 : Ns + Nrow, stab_size + 1 : Nrow];
        rows_scenario3 = tableau.tab(indices_scenario3, :);

        % Compute the symplectic inner product for all rows in scenario 3
        phases_scenario3 = Util.symplectic_inner_product_vec(rows_scenario3, row_measure, Ns);

        % Check if any row satisfies the condition
        idx_scenario3 = find(phases_scenario3 == 1, 1);
        if ~isempty(idx_scenario3)
            scenario = 3;
            row_idx = indices_scenario3(idx_scenario3);
            return;
        end
        
        % % scenario 1
        % for i = Ns+1: Ns + stab_size
        %     if Util.symplectic_inner_product(tableau.tab(i, :), row_measure, Ns)
        %         scenario = 1;
        %         row_idx = i;
        %         return
        %     end
        % end
        % 
        % % scenario 3
        % for i = [Ns + stab_size + 1: Ns + Nrow, stab_size + 1: Nrow]
        %     if Util.symplectic_inner_product(tableau.tab(i, :), row_measure, Ns)
        %         scenario = 3;
        %         row_idx = i;
        %         return
        %     end
        % end

        % scenario 2
        scenario = 2;
        row_idx = 0;
    end

end
end

