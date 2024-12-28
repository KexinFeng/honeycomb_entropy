function output = replica_average(varargin)
tictime = clock();
ip = inputParser;
ip.addParameter('T', 0.1);
ip.addParameter('cir', 32);
ip.addParameter('simDim', 30);
ip.addParameter('plotting', 0);
ip.addParameter('use_gpu', 0);
ip.addParameter('Nsample', 200);
ip.addParameter('rand_thrd', 1);
ip.addParameter('fun', @bipartite);
ip.addParameter('update', 1);
ip.addParameter('pick_vertex', 1);
ip.addParameter('kap', 0);
ip.addParameter('use_parfor', 1);
ip.addParameter('data_dir', []);


ip.parse(varargin{:});
pars = ip.Results;
global Nsample T cir kap
T = pars.T;
cir = pars.cir;
simDim = pars.simDim;
plotting = pars.plotting;
use_gpu = pars.use_gpu;
Nsample = pars.Nsample;
rand_thrd = pars.rand_thrd;
fun = pars.fun;
update = pars.update;
pick_vertex = pars.pick_vertex;
kap = pars.kap;
use_parfor = pars.use_parfor;
data_dir = pars.data_dir;

fprintf('T=%f\n', T);
%%
% cur_dir = pwd;
filepath = fileparts(mfilename('fullpath'));
% cd(filepath)
addpath(genpath([filepath, '/../../utils']));
addpath(genpath([filepath, '/../../utils2']));
addpath([filepath, '/../']);
addpath([filepath, '/library/']);
% cd(cur_dir);

fun_info = functions(fun);
fun_name = fun_info.function;

disp(fun_name)
fprintf('T=%f, cir=%d, kap=%.2f, simDim=%d, Nsample=%d\n',...
     T, cir, kap, simDim, Nsample);
for i = 1:length(varargin)
    disp(varargin{i});
end
%%

if strcmp(fun_name, 'attenuation')
    path = sprintf('./data_atten/cir_%d_kap_%dpct/', cir, floor(kap*100));
    name = [];
elseif strcmp(fun_name, 'spinRaman')
    path = sprintf('./data_%s/', fun_name);
    name = sprintf('%s', 'sR_');
else 
    path = sprintf('./data_%s/', fun_name);
    name = sprintf('%s_', fun_name);
end

if ~isempty(data_dir)
    path = [data_dir, path];
end
    
triexpf(~exist(path, 'dir'), {@mkdir, path}, {});
name = [name, sprintf('T_%.3f_cir_%d_kap_%.2f_simDim_%d_Nsample_%d', ...
         T, cir, kap, simDim, Nsample)];

                    
if ~update && exist([path, name, '.mat'], 'file')
    load([path, name, '.mat']);
    disp('file exists, loaded');
else 
    disp('update or file not exist');
    disp([path, name, '.mat'])
    output = {};
    if strcmp(fun_name, 'attenuation')
        [Ihat, Ierr, ~, nbar] = core_straSamp2(1, simDim^2, fun, ...
            'cir', cir,'T', T, 'num', simDim, 'cal_weight_DOS', 0,...
            'use_gpu', use_gpu, 'use_msi', 1, 'rand_thrd', rand_thrd, ...
            'pick_vertex', pick_vertex, 'kap', kap);
        Ihat = Ihat{1};
        Ierr = Ierr{1};
        alphas = reshape(Ihat, simDim, simDim);

    elseif strcmp(fun_name, 'spinRaman') || strcmp(fun_name, 'bubble') || strcmp(fun_name, 'rho')
        if strcmp(fun_name, 'spinRaman')
            num_cell = 3;
        elseif strcmp(fun_name, 'bubble')
            num_cell = 4*2;
        elseif strcmp(fun_name, 'rho')
            num_cell = 4;
        end
        [Ihat, Ierr, args, nbar] = core_straSamp2(num_cell, simDim, fun, ...
            'cir', cir,'T', T, 'num', simDim, 'kap', kap,...
            'use_gpu', use_gpu, 'rand_thrd', rand_thrd, 'plotting', 0, 'use_parfor', use_parfor);     
    end
    
    plot_tmp = plotting;
    clear plotting update
    save([path, name, '.mat']);
    fprintf('data saved:\n %s\n', [path, name]);
    fprintf('Time elapsed: %f s\n', etime(clock(),tictime));
    plotting = plot_tmp;
end

if ~exist('args', 'var')
    args = Omgs;
    plotting = 0;
end
output = {Ihat, Ierr, args, nbar};

%%
if plotting
    digits = regexp(name, '\d*pa', 'match');
    d = str2double(regexprep(digits, 'pa', ''));
    title_str = regexprep(name, {'_(?=\D)','_(?=\d)','\d*pct'}, {' ','=',num2str(d/100)});
    title_str = [title_str, sprintf(' Nphi = %.2g', nbar*cir*cir)];
    
    if strcmp(fun_name, 'attenuation')
        try
            title_str = [title_str, sprintf(' \nerr=%.3g', mean(Ierr))];
        catch
        end
        path = sprintf('./figures_atten/cir_%d_kap_%dpct/', cir, floor(kap*100));
        mkdir(path);
        save_str = [path, name];
        plot_alpha(alphas, title_str, save_str);
    else
        path = sprintf('./figures_%s/', fun_name);
        mkdir(path);
        save_str = [path, name];
        
        plot_intensity(Ihat, Omgs, title_str, save_str, 'err', Ierr, 'legend_str', {'all', 'pp', 'ph'});
    end
    fprintf('fig saved:\n %s\n', save_str);
end

fprintf('plot saved to %s\n', path);
fprintf('Time elapsed: %f s\n', etime(clock(),tictime));
end
 

%% embedding utils
function [Ihat, Ierr, Omgs, nbar] = core_straSamp2(num_cell, simDim, fun, varargin)
    global cir T Nsample kap
    
    wid = cir;
    info = functions(fun);
    fun_name = info.function;
    
    [probs, ~, ~,ns_raw, nbar, Efit] = prob_dist(cir, wid, T, 'kap', kap, 'plotting', 0);
    pause(0.001);

%     prob_cut = 1e-12/Nsample;
    prob_cut = 1e-2/Nsample;
    idxs = find(probs>=prob_cut);
    

    Ij_hat = cellfun(@(c) zeros(simDim, length(idxs)),...
        cell(1, num_cell), 'uni', false);
    Ij_var = cellfun(@(c) zeros(simDim, length(idxs)), ...
        cell(1, num_cell), 'uni', false);
    Omgs = [];


    nums = zeros(size(idxs));
    ps = probs(idxs)/sum(probs(idxs));
    ns = ns_raw(idxs);
    for j = 1:length(idxs)
%         nums(j) = round(Nsample * ps(j));
        min_sample = round(max(Nsample * 0.005, 1));
        nums(j) = max(round(Nsample * ps(j)), min_sample);
    end
    min_sample
    length(idxs)/length(probs)
    total_sample = sum(nums)
    fprintf('exceed budget number %d by %.2f\n', Nsample, (total_sample - Nsample)/Nsample);
    
    bws_bar = zeros(size(idxs));    
    for j = 1:length(idxs)
        n = ns(j);
        numj = nums(j);
        fprintf('accumulated samples=%.2f, percent = %.2f\n', sum(nums(1:j-1)), sum(nums(1:j-1))/total_sample);


        fbatch = cellfun(@(c) zeros(simDim, numj), ...
            cell(1, num_cell), 'uni', false);
        
        bws = zeros(1, numj);
        for s = 1:numj
            if strcmp(fun_name, 'attenuation')
                [yi, bwi] = fun(varargin{:}, 'nflux', n);
                Omgs = [];
            elseif strcmp(fun_name, 'spinRaman') || strcmp(fun_name, 'bubble')|| strcmp(fun_name, 'rho')
                [yi, bwi, Omgs] = fun(varargin{:}, 'nflux', n);
            end
            
            yi = triexp(~iscell(yi), {yi}, yi);
            
            try
                for c = 1:num_cell
                    fbatch{c}(:, s) = yi{c}(:); % all tensor linearized here
                end
            catch
                error('num_cell or simDim doesnt match fun output');
            end
            bws(s) = bwi;
        end
        fbar = cellfun(@(fb) fb.* reshape(bws, 1, [])/mean(bws), ...
            fbatch, 'uni', false);
        bws_bar(j) = mean(bws);
        
        for c = 1:num_cell
            Ij_hat{c}(:, j) = mean(fbar{c}, 2);
            Ij_var{c}(:, j) = std(fbar{c}, 0, 2).^2 / numj;
            if any(isnan(Ij_hat{c}), 'all')
                fprintf('cell = %d, j=%d\n', c, j);
                Ij_hat{c}
                bws
                disp(' ');
                error('NaN found in Ij_hat, possibly mean(bws) is zero');
            end
        end

%         sum(nums(1:j)) / Nsample
    end
    ps = ps .* bws_bar./exp(-Efit(idxs)/T); 
    % modify the model dist with a factor estimated with the samples
    
    ps = ps/sum(ps);
    if any(isnan(ps), 'all')
        ps
        bws_bar
        disp(' ');
        error('NaN found in ps');
    end
    
    Ihat = cellfun(@(Ij) sum(Ij .* reshape(ps, 1, []), 2), ...
        Ij_hat, 'uni', false);
    Ierr = cellfun(@(Ijv) sqrt(sum(Ijv .* reshape(ps.^2, 1, []), 2)),...
        Ij_var, 'uni', false);
end





