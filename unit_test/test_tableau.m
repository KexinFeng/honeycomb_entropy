filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
clc
rng(24);

%%
tab = load('./data/tab_c2w2T10.mat');
tableau = Tableau('cir', 2, 'wid', 2, ...
    'tab', tab.tab, 'stab_size', tab.stab_size);

g1 = Simulator();
g1.set_tableau(tableau);

g2 = Simulator();
g2.set_tableau(tableau);

tableau.render_table();
tableau.pair_tab_property();
g2.tableau.render_table();
g2.tableau.pair_tab_property();

%% Test partial trace
tableau = g1.tableau;
tableau.partial_trace(2);
assert(tableau.stab_size == 2);
tableau.render_table();
tableau.pair_tab_property();

tableau.partial_trace(5);
assert(tableau.stab_size == 1);
tableau.render_table();
tableau.pair_tab_property();

tableau2 = g2.tableau;
tableau2.partial_trace([2, 5]);
tableau1 = tableau;

assert(tableau2.stab_size == 1);
assert(all(tableau2.tab(1:tableau2.stab_size, :) == tableau1.tab(1:tableau1.stab_size, :), 'all'));
tableau2.render_table();
tableau2.pair_tab_property();


%%
tab = load('./data/tab_c2w2T30.mat');
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

%% Test partial trace
tableau = g1.tableau;
tableau.partial_trace(2);
assert(tableau.stab_size == 2);
tableau.render_table();
tableau.pair_tab_property();

tableau.partial_trace(5);
assert(tableau.stab_size == 1);
tableau.render_table();
tableau.pair_tab_property();

tableau2 = g2.tableau;
tableau2.partial_trace([2, 5]);
tableau1 = tableau;

% assert(tableau2.stab_size == 1);
assert(all(tableau2.tab(1:tableau2.stab_size, :) == tableau1.tab(1:tableau1.stab_size, :), 'all'));
tableau2.render_table();
tableau2.pair_tab_property();


%% Test clone
tableau_new = tableau1.clone();
tableau_new.tab(1) = 10;
assert(~all(abs(tableau_new.tab - tableau.tab)< 1e-3, 'all'))

%%
fprintf('success\n')
dbstop = 1;

