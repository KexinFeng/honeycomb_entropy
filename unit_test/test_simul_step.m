filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
rng(24);

%%
simul = Simulator('cir', 2, 'wid', 2);
cir = simul.cir;
wid = simul.wid;
Ns = cir * wid *2;

tableau = simul.tableau;
check_generator = simul.check_generator;

%% Y3Y5
row1 = check_generator.bond2row(1, 1, 3); % Y3Y5

[scenario, row_idx] = simul.check_scenario(row1);
assert(scenario == 3);

simul.scenario3(row1, row_idx);
assert(simul.tableau.stab_size == 1);


simul.tableau.render_table();
simul.tableau.pair_tab_property();


%% Z4Z7
row2 = check_generator.bond2row(2, 1, 1);

[scenario, row_idx] = simul.check_scenario(row2);
assert(scenario == 3);

simul.scenario3(row2, row_idx);
assert(tableau.stab_size == 2);

simul.tableau.render_table();
simul.tableau.pair_tab_property();

%% Y3Y5 * Z4Z7
row3 = Util.pauli_product(row1, row2);


[scenario, ~] = simul.check_scenario(row3);
assert(scenario == 2);

assert(tableau.stab_size == 2);

simul.tableau.render_table();
simul.tableau.pair_tab_property();

%% X3X7
row4 = check_generator.bond2row(2, 1, 2);

[scenario, row_idx] = simul.check_scenario(row4);
assert(scenario == 1);

simul.scenario1(row4, row_idx);
assert(tableau.stab_size == 2);


simul.tableau.render_table();
simul.tableau.pair_tab_property();

%% X2X6
row5 = check_generator.bond2row(1, 2, 2);

[scenario, row_idx] = simul.check_scenario(row5);
assert(scenario == 3);

simul.scenario3(row5, row_idx);
assert(tableau.stab_size == 3);


simul.tableau.render_table();
simul.tableau.pair_tab_property();

%%
fprintf('success\n')
dbstop = 1;