function [] = measure_simul(varargin)
    tic 
    filepath = fileparts(mfilename('fullpath'));
    addpath(genpath([filepath, '/../utils']));
    % addpath([filepath, '/../']);
    
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
            dbstop = 1;
        end
        fprintf('step=%d, r=%d\n', step, stab_size)

        render_table(tableau, stab_size, cir, wid);
        if ~pair_commute(tableau, stab_size, cir, wid)
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


function [row, bond] = generate_bond(x, y, pars)
    cir = pars.cir;
    wid = pars.wid;
    Ns = cir*wid*2;
    boundary = pars.boundary;

    row = zeros(1, Ns*2);

    % Site index: x = 1~cir; y = 1~wid; ab = 0,1
    pos = @(x, y, ab) y + wid.*(x-1) + cir*wid.*ab;
    
    % Generate bond
    if strcmp(boundary, "periodic")
        bond = randi([1, 3]); % 1:z, 2:x, 3:y
    elseif strcmp(boundary, 'open')
        if x==1 && y==1 
            bond = randsample([1, 3], 1);
        elseif x == cir && y == wid
            bond = [];
            return
        elseif y == wid
            bond = randsample([2, 3], 1);
        elseif x == cir
            bond = randsample([1, 2], 1);
        else
            bond = randi([1, 3]);
        end
    else
        error('%s not right', boundary)
    end
    
    % (x, y, 0)-X-(x, y , 1)
    % (x, y+1, 0)-Z-(x, y, 1)
    % (x+1, y, 0)-Y-(x, y, 1)
    bin = dec2bin(bond, 2);
    bin_arr = kron(str2num(bin(:)), [1; 1]);
    
    xs = [x, homod(x+(bond == 3), cir)];
    ys = [y, homod(y+(bond == 1), wid)];
    ab = [1, 0];
    
    idx = pos(xs, ys, ab);
    idx = [idx, idx + Ns];
    
    row(idx) = reshape(bin_arr, size(idx));
    % bond_names = ['Z', 'X', 'Y'];
    % pauli = row2pauli(row, cir, wid);
    % fprintf('%s\n pauli: %s\n', bond_names{bond}, pauli)
end


function [tab, stab_size, bond] = measure_rank(x, y, tab, stab_size, pars)
    cir = pars.cir;
    wid = pars.wid;
    Ns = cir*wid*2;

    [row, bond] = generate_bond(x, y, pars);
    
    % Find the first anti-commuting row
    phase = 0;
    idx = 1;
    while idx <= stab_size
        phase = symplectic_inner_product(row, tab(idx, :), Ns);
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
           phase = symplectic_inner_product(row, tab(idx, :), Ns);
           if phase == 1
               tab(idx, :) = pauli_product(row_discard, tab(idx, :));
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

%% utility functions
function phase = symplectic_inner_product(row1, row2, Ns)
    assert(size(row1, 2) == size(row2, 2) && size(row2, 2) == Ns*2);
    phase = sum(row1(1, 1:Ns).* row2(1, Ns+1:end) + row2(1, 1:Ns).* row1(1, Ns+1:end));
    phase = mod(phase, 2);
end

function row = pauli_product(row1, row2)
    assert(all(size(row1) == size(row2)));
    row = mod(row1 + row2, 2);
end

function [] = render_table(tab, r, cir, wid)
    for i = 1:r
        fprintf("%s\n", row2pauli(tab(i, :), cir, wid));
    end
end

function b = pair_commute(tab, r, cir, wid)
    b = true;
    for i = 1:r
        for j = i+1:r
            b = 0 == symplectic_inner_product(tab(i, :), tab(j, :), cir*wid*2);
            if ~b
                break
            end
        end
    end
end

function pauli_str = row2pauli(row, cir, wid)
    % 2:x, 1:z, 3:y
    bond_names = ['Z', 'X', 'Y'];
    Ns = cir * wid *2;
    strArray = [];
    assert(size(row, 2) == Ns * 2);
    for ind = 1: Ns
        bin = sum(row([ind, ind+Ns]).*[2, 1]);
        if bin == 0
            continue;
        end
        operator = bond_names(bin);
        % pos = @(x, y, ab) y + wid.*(x-1) + cir*wid.*ab;
        [y, x, a] = ind2sub([wid, cir, 2], ind);
        % [res, a1] = homod(ind, cir*wid);
        % [res, x1] = homod(res, wid);
        % y1 = res;
        % assert(x==x1 && y==y1 && a==a1);
        a = a-1;
        ab = triexp(a==0, "a", "b");
        string = sprintf("%s_%d", operator, ind);
        strArray = [strArray, string];
    end
    pauli_str = triexp(~isempty(strArray), strjoin(strArray, " "), "");
end

function [] = nothing()
%%
vertx_str = {'MS', 'E2', 'A1'};
output = cell(1, 3);
M = triexp(use_msi, maxNumCompThreads(), 1);

% rng(1);
rng('shuffle');
numplaq = cir * cir;

fluxes = unique(fluxes);
if isempty(fluxes)
    fluxes = randperm(numplaq, nflux);
else
    nflux = length(fluxes);
    if mod(nflux, 2) ~= 0
        warning('odd number of fluxes will be trimed');
        nflux = nflux -1;
        fluxes(end) = [];
    end
end
fluxes = sort(fluxes);
thrd_state = triexp(rand_thrd, randperm(4, 1)-2, thrd_state);

% delta = .1;
qxlim = .5*pi;
qylim = .5*pi;
qx_coor = linspace(-qxlim, qxlim, num);
qy_coor = linspace(-qylim, qylim, num);
qx_coor = qx_coor + 1e-9 * rand(1, num);
qy_coor = qy_coor + 1e-9 * rand(1, num);

% assert_qnorm(qx_coor, qy_coor)
% assert_thetaq(qx_coor, qy_coor, mu)
c1 = 1.5;
c2 = 1.5;
rho = 1;
assert_vs_vF(c1, c2, rho);

wid = cir;
Nx = cir;
N = cir * wid;


H = hamiltonianT(cir, wid, 0, 'fluxes', fluxes, 'change_gauge', thrd_state,...
    'plotting', plotting, 'kap', kap);
if plotting
    path = sprintf('./figures_app/%d/', cir);
    mkdir(path);
    save_str = [path, sprintf('%d_%d_nflux_%d_fluxes_%s', cir, wid, nflux,  DataHash(fluxes))];
    print(gcf, '-depsc2', [save_str,'.eps']);
    eps2pdf([save_str,'.eps'],[save_str,'.pdf'],1);
    delete([save_str, '.eps']);
end
[u, s] = diagnlz2(H, [cir, wid], kap, 'use_gpu', use_gpu);
s = double(gather(s));
u = double(gather(u));
bw = boltzmann_weight_cal(T, s, cir, wid, 'kap', kap);

alphas = zeros(num, num);
simSpace = [num, num];

% % This potentially prevents the intialization problem in parpool
% launch_parpool();
% parfor idx = 1:prod(simSpace)
for idx = 1:prod(simSpace)
    [i, j] = ind2sub(simSpace, idx);    %  [i, j] = ind2sub([3, 7], 1:21);
    qx = qx_coor(i);
    qy = qy_coor(j);
    qxy = [qx, qy];
    qxy = qxyclip(qxy);
    qxy = nearest_onsite(qxy, [cir, wid]);

%     c1 = .01;
%     c2 = .01;
%     rho = 1;
    [rq, vs] = phonon_mode(c1, c2, rho, qxy);
    vs(1) = 1.715848;

    Omgs = vs * norm(qxy);
    if strcmp(pick_core, 'TS')
        [Piph, Pipp, ~] = bubbleTS2(qxy, Omgs(1), delta, Nx, T, ...
            'onshell', 0, 'mu', mu, 'pick_vertex', pick_vertex);
    else
        [Piph, Pipp, ~] = bubbleDS2(qxy, Omgs(1), delta, cir, T,...
            'onshell', 0, 'mu', mu, 'pick_vertex', pick_vertex, 'use_gpu', use_gpu,...
            'u', u, 's', s, 'fluxes', fluxes, 'thrd_state', thrd_state);
    end
    Bphs = imag(Piph);
    Bpps = imag(Pipp);

    if strcmp(channel, 'ph')
        Omg_dist = Bphs;
    elseif strcmp(channel, 'phpp')
        Omg_dist = Bphs + Bpps;
    end

    % alpha = -1/(vs^2 * q) imag(Pi(q, Omg=vs*q))
    alpha = triexp(denoise(Omg_dist, 1e-20)~=0, 1./vs(mu).^2 ./ norm(qxy) .* Omg_dist, 0);
    alphas(idx) = double(gather(alpha));
end
output{1} =  inversion_avg(alphas);
toc



%%
if plotting
    title_str = sprintf('%s channel = %s, delta = %.2f, T = %.2f, N_x = %d,mu=%d, gauge=%d\n num=%d qlim=[%.2f, %.2f]\\pi',...
        pick_core, channel, delta, T, cir, mu, thrd_state, num, qxlim/pi, qylim/pi);
    path = sprintf('./figures_app/%d_T_%dpct/', cir, floor(T*100));
    mkdir(path);
    save_str = [path, sprintf('%s', triexp(use_msi, 'msi_', ''))];
    save_str = [save_str, sprintf('%s_channel_%s_T_%dpa_cir_%d_mu_%d_thrd_%d_qlim_%dpa_%dpa_num_%d_delta_%dpa', ...
        pick_core, channel, floor(T*100), cir, mu, thrd_state, ...
        floor(qxlim/pi*100), floor(qylim/pi*100), num, floor(delta*100))];
    save_str = [save_str, triexp(nflux>0, sprintf('_nflux_%d_fluxes_%s', nflux,  DataHash(fluxes)), '')];
    title_str = [title_str, triexp(nflux>0, sprintf(' nflux=%d', nflux), '')];
    plot_alpha(alphas, title_str, save_str);

    mkdir('./data_app/');
    save(['./data_app/',strrep(save_str, path , '')]);
end
end



