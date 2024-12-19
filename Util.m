classdef Util
    methods(Static)
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
                fprintf("%s\n", Util.row2pauli(tab(i, :), cir, wid));
            end
        end

        function [] = render_table_destab(tab, stab_size, cir, wid)
            Ns = cir * wid * 2;
            for i = 1 : Ns*2
            % for i = 2 : Ns*2-1
            %     if i == Ns || i == Ns + 1
            %         continue
            %     end
                fprintf("%d: %s\n", i, Util.row2pauli(tab(i, :), cir, wid));
            end
        end

        function b = pair_tab_property(tab, r, cir, wid)
            Ns = cir * wid * 2;
            b = true;
            for i = 1:Ns*2
                for j = i+1:Ns*2
                    if j - Ns == i
                        b = 1 == Util.symplectic_inner_product(tab(i, :), tab(j, :), Ns);
                    else
                        b = 0 == Util.symplectic_inner_product(tab(i, :), tab(j, :), Ns);
                    end
                    if ~b
                        break
                    end
                end
            end
        end
        
        function b = pair_commute(tab, r, cir, wid)
            b = true;
            for i = 1:r
                for j = i+1:r
                    b = 0 == Util.symplectic_inner_product(tab(i, :), tab(j, :), cir*wid*2);
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
            if ~isempty(strArray)
                pauli_str = strjoin(strArray, " ");
            else
                pauli_str = "";
            end
        end

    end
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

