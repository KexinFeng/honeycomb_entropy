filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));

simul = SimulatorDestab('cir', 2, 'wid', 2, 'boundary', 'open');
simul.verbose = true;

cir = simul.cir;
wid = simul.wid;
Ns = cir * wid *2;

tab = eye(Ns*2, Ns*2);
stab_size = 0; % stablizer size
if strcmp(simul.boundary, 'open')
    tab([1, Ns+1, Ns-1, 2*Ns-1], :) = tab([Ns-1, 2*Ns-1, 1, Ns+1], :);
    % print
    Util.render_table_destab(tab, 0, cir, wid);
end

%% Y3Y5
row1 = simul.bond2row(1, 1, 3); % Y3Y5

[scenario, row_idx] = simul.check_scenario(tab, row1, stab_size);
assert(scenario == 3);

[tab, stab_size] = simul.scenario3(tab, row1, row_idx, stab_size);
assert(stab_size == 1);

Util.render_table_destab(tab, stab_size, cir, wid);
Util.pair_tab_property(tab, cir, wid);

dbstop = 1;

%% Z4Z6
row2 = simul.bond2row(2, 1, 1);

[scenario, row_idx] = simul.check_scenario(tab, row2, stab_size);
assert(scenario == 3);

[tab, stab_size] = simul.scenario3(tab, row2, row_idx, stab_size);
assert(stab_size == 2);

Util.render_table_destab(tab, stab_size, cir, wid);
Util.pair_tab_property(tab, cir, wid);

%% Y3Y5 * Z4Z6
row3 = Util.pauli_product(row1, row2);
% % print
% Util.row2pauli(row3, cir, wid)

[scenario, ~] = simul.check_scenario(tab, row3, stab_size);
assert(scenario == 2);

assert(stab_size == 2);

Util.render_table_destab(tab, stab_size, cir, wid);
Util.pair_tab_property(tab, cir, wid);

%% X3X6
row4 = simul.bond2row(2, 1, 2);

[scenario, row_idx] = simul.check_scenario(tab, row4, stab_size);
assert(scenario == 1);

tab = simul.scenario1(tab, row4, row_idx);
assert(stab_size == 2);

Util.render_table_destab(tab, stab_size, cir, wid);
Util.pair_tab_property(tab, cir, wid);

%%
fprintf('\nsuccess\n')
dbstop = 1;