function [] = measure_simul_rank(varargin)
    tic 
    filepath = fileparts(mfilename('fullpath'));
    addpath(genpath([filepath, '/../utils']));
    % addpath([filepath, '/util.m']);
    
    rng(24);
    
    %% period boundary
    ip = inputParser;
    ip.addParameter('cir', 2);
    ip.addParameter('wid', 2);
    ip.addParameter('plotting', 1);
    ip.addParameter('boundary', 'open');
    ip.addParameter('T', 10);
    
    ip.parse(varargin{:});
    pars = ip.Results;
    
    cir = pars.cir;
    wid = pars.wid;
    plotting = pars.plotting;
    boundary = pars.boundary;
    
    %%
    Ns = cir*wid*2;
    
    tableau = zeros(Ns+1, Ns*2);
    stab_size = 0; % stablizer size
    
    T = pars.T;
    
    if plotting
        figure
        hold on
        plot(0, Ns - 2*strcmp(boundary, 'open') - stab_size, '*', 'Color', 'b');
        pause(0.01);        
    end
    
    for step = 1:T
        for idx = 1:cir*wid
            [y, x] = homod(idx, wid);
            [tableau, stab_size, bond] = measure_rank(x, y, tableau, stab_size, pars);
            % [tableau, stab_size, bond] = measure_destab(x, y, tableau, stab_size, pars);
            dbstop = 1;
        end

        % Printing
        fprintf('step=%d, r=%d\n', step, stab_size)
        Util.render_table(tableau, stab_size, cir, wid);
        if ~Util.pair_commute(tableau, stab_size, cir, wid)
            error('not commute')
        end
        disp(' ')

        % Plotting
        if plotting
            plot(step, Ns - 2*strcmp(boundary, 'open') - stab_size, '*', 'Color', 'b');
            pause(0.01);
        end
    end
    
    disp(['r=', num2str(stab_size), ' total spin:', num2str(Ns - 2*strcmp(boundary, 'open'))])

end


function [tab, stab_size, bond] = measure_rank(x, y, tab, stab_size, pars)
    cir = pars.cir;
    wid = pars.wid;
    Ns = cir*wid*2;

    [row, bond] = generate_bond(x, y, pars);

    if isempty(bond)
        return
    end
    
    % Find the first anti-commuting row
    phase = 0;
    idx = 1;
    while idx <= stab_size
        phase = Util.symplectic_inner_product(row, tab(idx, :), Ns);
        if phase  == 1
            break
        end
        idx = idx + 1;
    end
    
    if phase == 1
        % Swap it with the new Pauli
        row_discard = tab(idx, :);
        tab(idx, :) = row;
        idx = idx + 1;
        % Modify the rest anti-commuting rows
        while idx <= stab_size
           phase = Util.symplectic_inner_product(row, tab(idx, :), Ns);
           if phase == 1
               tab(idx, :) = Util.pauli_product(row_discard, tab(idx, :));
           end
           idx = idx + 1;
        end
    else
        % New stablizer is added if row is independent from existing ones
        tab(stab_size+1, :) = row;
        tab_galois = gf(tab(1: stab_size+1, :), 1, 2);
        if rank(tab_galois) == stab_size + 1
            stab_size = stab_size + 1;
        end
    end
end




