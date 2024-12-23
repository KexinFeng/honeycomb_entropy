filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
clc

%%
tab = load('/Users/kx/Desktop/forked/K_circuit/code/honeycomb_entropy/unit_test/tab_c2w2T10.mat');
tableau = Tableau('cir', 2, 'wid', 2, ...
    'tab', tab.tab, 'stab_size', tab.stab_size);

g1 = Simulator();
g1.set_tableau(tableau);

g2 = Simulator();
g2.set_tableau(tableau);

g1.tableau.render_table();
g1.tableau.pair_tab_property();
g2.tableau.render_table();
g2.tableau.pair_tab_property();

%%
g1.partial_trace(2);
assert(g1.tableau.stab_size == 2);
g1.tableau.render_table();
g1.tableau.pair_tab_property();

g1.partial_trace(5);
assert(g1.tableau.stab_size == 1);
g1.tableau.render_table();
g1.tableau.pair_tab_property();

g2.partial_trace([2, 5]);
tableau2 = g2.tableau;
tableau1 = g1.tableau;

assert(tableau2.stab_size == 1);
assert(all(tableau2.tab(1:tableau2.stab_size, :) == tableau1.tab(1:tableau1.stab_size, :), 'all'));
g2.tableau.render_table();
g2.tableau.pair_tab_property();


%%
tab2 = load('/Users/kx/Desktop/forked/K_circuit/code/honeycomb_entropy/unit_test/tab_c2w2T30.mat');
tableau = Tableau('cir', 2, 'wid', 2, ...
    'tab', tab.tab, 'stab_size', tab.stab_size);

g1 = Simulator();
g1.set_tableau(tableau);

g2 = Simulator();
g2.set_tableau(tableau);

g1.tableau.render_table();
g1.tableau.pair_tab_property();
g2.tableau.render_table();
g2.tableau.pair_tab_property();

%%
g1.partial_trace(2);
assert(g1.tableau.stab_size == 2);
g1.tableau.render_table();
g1.tableau.pair_tab_property();

g1.partial_trace(5);
assert(g1.tableau.stab_size == 1);
g1.tableau.render_table();
g1.tableau.pair_tab_property();

g2.partial_trace([2, 5]);
tableau2 = g2.tableau;
tableau1 = g1.tableau;

% assert(tableau2.stab_size == 1);
assert(all(tableau2.tab(1:tableau2.stab_size, :) == tableau1.tab(1:tableau1.stab_size, :), 'all'));
g2.tableau.render_table();
g2.tableau.pair_tab_property();


%%
fprintf('success\n')
dbstop = 1;
