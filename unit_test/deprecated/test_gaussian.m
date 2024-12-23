filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
clc

tab = load('../tmp/tab4unit_test.mat');
tableau = tab.tableau;
stab_size = tab.stab_size;


g1 = GaussianEliminator('cir', 2, 'wid', 2, ...
    'tab', tableau, 'stab_size', stab_size);

g2 = GaussianEliminator('cir', 2, 'wid', 2, ...
    'tab', tableau, 'stab_size', stab_size);

%%
Util.render_table_destab(g1.tab, g1.stab_size, 2, 2);

g1.partial_trace(2);
assert(g1.stab_size == 2);
Util.render_table_destab(g1.tab, g1.stab_size, 2, 2);
Util.pair_tab_property(g1.tab, 2, 2)

g1.partial_trace(5);
assert(g1.stab_size == 1);
Util.render_table_destab(g1.tab, g1.stab_size, 2, 2);

g2.partial_trace([2, 5]);
assert(g2.stab_size == 1);
assert(all(g2.tab(1:g2.stab_size, :) == g1.tab(1:g1.stab_size, :), 'all'));
Util.render_table_destab(g2.tab, g2.stab_size, 2, 2);

Util.pair_tab_property(g1.tab, 2, 2);

% dbstop = 1;
