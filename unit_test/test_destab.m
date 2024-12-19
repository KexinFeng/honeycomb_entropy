%%
filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
rng(24);
clc;

%%
pars.cir = 2;
pars.wid = 2;
pars.boundary = 'open';

cir = pars.cir;
wid = pars.wid;
Ns = pars.cir * pars.wid *2;

tab = eye(Ns*2, Ns*2);
stab_size = 0; % stablizer size
if strcmp(pars.boundary, 'open')
    tab([1, Ns+1, Ns-1, 2*Ns-1], :) = tab([Ns-1, 2*Ns-1, 1, Ns+1], :);
    % print
    Util.render_table_destab(tab, 0, cir, wid);
end

%% get to-test func
input_cell = {'testing', true, 'plotting', false, 'cir', cir, 'wid', wid};
res = measure_simul_destab(input_cell{:});

for func = res.test_funcs
    name = func2str(func{1});
    fprintf('%s\n', name);
end

check_scenario = res.test_funcs{3};
scenario1 = res.test_funcs{1};
scenario3 = res.test_funcs{2};
[~, ~, res] = generate_bond(1, 1, pars);

%% Y3Y5
row1 = res.bond2row(1, 1, 3, pars); % Y3Y5

[scenario, row_idx] = check_scenario(tab, row1, stab_size, pars);
assert(scenario == 3);

[tab, stab_size] = scenario3(tab, row1, row_idx, stab_size);
assert(stab_size == 1);

Util.pair_tab_property(tab, cir, wid);

%% Z4Z6
row2 = res.bond2row(2, 1, 1, pars);

[scenario, row_idx] = check_scenario(tab, row2, stab_size, pars);
assert(scenario == 3);

[tab, stab_size] = scenario3(tab, row2, row_idx, stab_size);
assert(stab_size == 2);

Util.render_table_destab(tab, stab_size, cir, wid);
Util.pair_tab_property(tab, cir, wid);

%% Y3Y5 * Z4Z6
row3 = Util.pauli_product(row1, row2);
Util.row2pauli(row3, cir, wid)

[scenario, ~] = check_scenario(tab, row3, stab_size, pars);
assert(scenario == 2);

Util.pair_tab_property(tab, cir, wid);
%%
Util.render_table_destab(tab, stab_size, cir, wid);


dbstop = 1
    
