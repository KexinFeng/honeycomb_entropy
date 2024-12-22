filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
clc

tab = load('../tmp/tab4unit_test.mat');
tableau = Tableau('cir', 2, 'wid', 2, ...
    'tab', tab.tableau, 'stab_size', tab.stab_size);

g1 = Simulator();
g1.set_tableau(tableau);

g2 = Simulator();
g2.set_tableau(tableau);

%%
% Util.render_table_destab(g1.tab, g1.stab_size, 2, 2);
g1.tableau.render_table_destab();

g1.partial_trace(2);
assert(g1.tableau.stab_size == 2);
g1.tableau.render_table_destab();
g1.tableau.pair_tab_property();
% Util.render_table_destab(g1.tab, g1.stab_size, 2, 2);
% Util.pair_tab_property(g1.tab, 2, 2)

g1.partial_trace(5);
assert(g1.tableau.stab_size == 1);
% Util.render_table_destab(g1.tab, g1.stab_size, 2, 2);
g1.tableau.render_table_destab();
g1.tableau.pair_tab_property();

g2.partial_trace([2, 5]);
tableau2 = g2.tableau;
tableau1 = g1.tableau;

assert(tableau2.stab_size == 1);
assert(all(tableau2.tab(1:tableau2.stab_size, :) == tableau1.tab(1:tableau1.stab_size, :), 'all'));
g2.tableau.render_table_destab();
g2.tableau.pair_tab_property();


% dbstop = 1;
