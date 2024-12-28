filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
rng(24);

%% identifier
idnt = Identifier('cir', 7, 'wid', 6, 'shift', 'rect');
wid = idnt.wid;

% single check
[x, y] = idnt.comb(4, 0);
assert(x == 1)
assert(y == 6)

delta_ys_a_1 = 0: 2: wid - 2;
delta_xs_a_1 = 0: -1: -(wid/2 - 1);
delta_ys_b_1 = 1: 2: wid - 1;
delta_xs_b_1 = -1: -1: -wid/2;

ys = 1 + [delta_ys_a_1, delta_ys_b_1];
xs = 1 + [delta_xs_a_1, delta_xs_b_1];
as = [zeros(size(delta_ys_a_1)), ones(size(delta_ys_b_1))];

[xs_out, ys_out] = idnt.comb(xs, ys);

output_exp_0 = [...
1 1 0
7 3 0
6 5 0
7 2 1
6 4 1
5 6 1];

for i = 1: length(xs_out)
    fprintf('%d %d %d\n', xs_out(i), ys_out(i), as(i))
    assert(all([xs_out(i), ys_out(i), as(i)] == output_exp_0(i, :), 'all'));
end
fprintf('\n');

%% benchmark len2qubits
check_generator = CheckGenerator('cir', 7, 'wid', 6);

linidx1 = check_generator.len2qubits(3);
linidx2 = check_generator.len2qubits(3, 'start', 2);
lindex = setdiff(linidx1, linidx2);
[xs, ys, as] = check_generator.pos_reverse(lindex);

output_exp = [...
1 1 1
1 2 1
6 5 1
6 6 1
7 3 1
7 4 1
1 1 2
5 6 2
6 4 2
6 5 2
7 2 2
7 3 2];

for i = 1: length(xs)
    fprintf('%d %d %d\n', xs(i), ys(i), as(i))
    assert(all([xs(i), ys(i), as(i)] == output_exp(i, :), 'all'));
end

fprintf('\n');

%% benchmark len2qubits with rect shift
check_generator = CheckGenerator('cir', 7, 'wid', 6, 'shift', 'rect');

linidx1 = check_generator.len2qubits(3);
linidx2 = check_generator.len2qubits(3, 'start', 2);
lindex2 = setdiff(linidx1, linidx2);
[xs, ys, as] = check_generator.pos_reverse(lindex2);

for i = 1: length(xs)
    fprintf('%d %d %d\n', xs(i), ys(i), as(i))
    assert(all([xs(i), ys(i), as(i)] == output_exp(i, :), 'all'));
end

fprintf('\n');

%%
fprintf('success\n')
dbstop = 1;

