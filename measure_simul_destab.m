function res = measure_simul_destab(varargin)
    tic 
    filepath = fileparts(mfilename('fullpath'));
    addpath(genpath([filepath, '/../utils']));
    % addpath([filepath, '/util.m']);
    
    rng(24);
    clc;
    res = struct();

    %% period boundary
    ip = inputParser;
    ip.addParameter('cir', 2);
    ip.addParameter('wid', 2);
    ip.addParameter('plotting', 1);
    ip.addParameter('boundary', 'open');
    ip.addParameter('T', 20);
    ip.addParameter('testing', false);
    
    ip.parse(varargin{:});
    pars = ip.Results;
    
    cir = pars.cir;
    wid = pars.wid;
    boundary = pars.boundary;
    
    %% prepare
    Ns = cir*wid*2;
    
    tableau = eye(Ns*2, Ns*2);
    if strcmp(pars.boundary, 'open')
        tableau([1, Ns+1, Ns-1, 2*Ns-1], :) = tableau([Ns-1, 2*Ns-1, 1, Ns+1], :);
        % print
        disp('swapped')
    end
    stab_size = 0; % stablizer size
    Util.render_table_destab(tableau, stab_size, cir, wid);

    T = pars.T;
    
    if pars.plotting
        figure
        hold on
        plot(0, Ns - 2*strcmp(boundary, 'open') - stab_size, '*', 'Color', 'b');
        pause(0.01);        
    end
    
    if pars.testing
        [~, ~, ~, res] = measure_destab(1, 1, tableau, stab_size, pars);
        return
    end

    %% main loop
    for step = 1:T
        for idx = 1:cir*wid
            [y, x] = homod(idx, wid);
            % [tableau, stab_size, bond] = measure_rank(x, y, tableau, stab_size, pars);
            [tableau, stab_size, bond, res] = measure_destab(x, y, tableau, stab_size, pars);
            if isempty(bond)
                continue
            end

            if ~Util.pair_tab_property(tableau, cir, wid)
                error('not commute')
            end
            % print
            fprintf('\nidx: %d\n', idx);
            Util.render_table_destab(tableau, stab_size, cir, wid);
        end

        if ~Util.pair_tab_property(tableau, cir, wid)
            error('not commute')
        end
        % print
        fprintf('step=%d, r=%d\n\n', step, stab_size)
        disp(' ')

        % Plotting
        if pars.plotting
            plot(step, Ns - 2*strcmp(boundary, 'open') - stab_size, '*', 'Color', 'b');
            pause(0.01);
        end
    end
    
    disp(['stab_size=', num2str(stab_size), ' total spin:', num2str(Ns - 2*strcmp(boundary, 'open'))])

end


function [tab, stab_size, bond, res] = measure_destab(x, y, tab, stab_size, pars)
    % Measure on the unit cell (x, y) 
    res = struct();

    cir = pars.cir;
    wid = pars.wid;
    Ns = cir*wid*2;


    [row, bond] = generate_bond(x, y, pars);
    % print
    if isempty(bond)
        return
    end
    fprintf('bond: %s\n', Util.row2pauli(row, cir, wid))

    [scenario, row_idx] = check_scenario(tab, row, stab_size, pars);
    
    if scenario == 1
        % stochastic output containing de-stablizer
        tab = scenario1(tab, row, row_idx);
    elseif scenario == 2
        % deterministic output
        % ;
    else
        % stochastic output containing enhanced space
        [tab, stab_size] = scenario3(tab, row, row_idx, stab_size);
    end

    if pars.testing
        res.test_funcs = {@scenario1, @scenario3, @check_scenario};
    end      
end


%% Functions
function tab = scenario1(tab, row_measure, row_idx)
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

function [tab, stab_size] = scenario3(tab, row_measure, row_idx, stab_size)
    % row_idx points to the row anticommuting with row_measure
    Ns = size(row_measure, 2) / 2;

    % swap row_idx to Ns + stab_size + 1
    row_idx_bar = homod(row_idx + Ns, 2*Ns);
    tab([row_idx, Ns + stab_size + 1], :) = tab([Ns + stab_size + 1, row_idx], :);
    tab([row_idx_bar, stab_size + 1], :) = tab([stab_size + 1, row_idx_bar], :);
    row_idx = Ns + stab_size + 1;

    tab = scenario1(tab, row_measure, row_idx);
    stab_size = stab_size + 1;
end


function [scenario, row_idx] = check_scenario(tab, row_measure, stab_size, pars)
    Ns = pars.cir * pars.wid * 2;
    
    Nrow = Ns - triexp(strcmp(pars.boundary, 'open'), 2, 0);

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


