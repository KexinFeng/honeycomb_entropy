function [output, xs] = bipartite(varargin)
tictime = clock(); 
script_path = fileparts(mfilename('fullpath'));
addpath(genpath([script_path, '/../../utils']));
addpath([script_path, '/../']);


%% input
ip = inputParser;
ip.KeepUnmatched = true;
ip.PartialMatching = false;

ip.addParameter('cir', 10);
ip.addParameter('wid', 8);
ip.addParameter('boundary', 'periodic');
ip.addParameter('shift', 0);
ip.addParameter('T', 150);
ip.addParameter('verbose', false);
ip.addParameter('probs', [0.25, 0.25, 0.25, 0.25]);
ip.addParameter('num', 10);
ip.addParameter('plotting', 1);
ip.addParameter('update', 0);

ip.parse(varargin{:});
pars = ip.Results;

pars.wid = pars.cir;
% Convert pars to a cell array of parameter-value pairs
fields = fieldnames(pars);      
values = struct2cell(pars);     
varargin = reshape([fields, values]', 1, []); 

%%
folder = sprintf([script_path, '/data_%s/'], 'bipartite');
triexpf(~exist(folder, 'dir'), {@mkdir, folder}, {});
name = sprintf('cir_%d_T_%d_probs_%.2f_%.2f_%.2f_%.2f_%s_%s',...
    pars.cir, pars.T, pars.probs(1), pars.probs(2), pars.probs(3), ...
    pars.probs(4), pars.boundary, pars.shift);

fprintf('%s\n', name);                
if ~ pars.update && exist([folder, name, '.mat'], 'file')
    load([folder, name, '.mat']);
    disp('file exists, loaded');
else
    disp('update or file not exist');
    
    %%
    simul = Simulator(varargin{:}, 'plotting', 0, 'verbose', false);
    [es, ts] = simul.simulate();
    
    %% save
    pars_tmp = pars;
    clear pars
    save([folder, name, '.mat']);
    fprintf('data saved:\n %s\n', [folder, name]);
    fprintf('Time elapsed: %f s\n',  etime(clock(), tictime));
    pars = pars_tmp;
end


%% plotting
if pars.plotting
    % title_str = {'title', 'sub'};    
    title_str = name;
    title_str = regexprep(title_str, '(?<=\D)_', '=');
    title_str = regexprep(title_str, '(?<=\d)_', ' ');
    
    folder = [script_path, sprintf('/figures_%s/', 'bipartite')];
    mkdir(folder);
    save_str = [folder, 'purify_', name];  
    plot_purify(es, ts, title_str, save_str, ...
        'legend_str', {sprintf('L=%d', pars.cir)});
end

%% measure
sys_sizes = 0: ceil(pars.cir / pars.num): pars.cir;
entropies = zeros(size(sys_sizes));
for idx = 1: length(sys_sizes)
    tableau = simul.tableau.clone();
    sys_size = sys_sizes(idx);
    qubits_env = simul.check_generator.len2qubits(pars.cir - sys_size, 'start', 1 + sys_size);
    tableau.partial_trace(qubits_env);
    entropies(idx) = tableau.get_entropy();
end

output = {entropies / pars.cir};
xs = sys_sizes / pars.cir;

%% plotting
if pars.plotting    
    % title_str = {'title', 'sub'};    
    title_str = name;
    title_str = regexprep(title_str, '(?<=\D)_', '=');
    title_str = regexprep(title_str, '(?<=\d)_', ' ');
    
    folder = [script_path, '/figures_bipartite/'];
    mkdir(folder);
    save_str = [folder, '/bipartite_', name];  
    plot_bipartite(output, xs, title_str, save_str, ...
        'legend_str', {sprintf('L=%d', pars.cir)});
end

end
