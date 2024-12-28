function [output, xs] = bipartite(varargin)
tictime = clock(); 
filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../../utils']));
addpath([filepath, '/../']);

rng(24);
clc;

%% input
ip = inputParser;
ip.KeepUnmatched = true;
ip.PartialMatching = false;

ip.addParameter('cir', 8);
ip.addParameter('wid', 8);
ip.addParameter('boundary', 'periodic');
ip.addParameter('shift', 0);
ip.addParameter('T', 30);
ip.addParameter('verbose', false);
ip.addParameter('probs', [1/4, 1/4, 1/4, 1/4]);
ip.addParameter('num', 8);
ip.addParameter('plotting', 1);
ip.addParameter('update', 1);

ip.parse(varargin{:});
pars = ip.Results;
pars.wid = pars.cir;

% Convert pars to a cell array of parameter-value pairs
fields = fieldnames(pars);      
values = struct2cell(pars);     
varargin = reshape([fields, values]', 1, []); 

%%
path = sprintf('./data_%s/', 'bipartite');
triexpf(~exist(path, 'dir'), {@mkdir, path}, {});
name = sprintf('cir_%d_T_%d_probs_%.2f_%.2f_%.2f_%.2f_%s',...
    pars.cir, pars.T, pars.probs(1), pars.probs(2), pars.probs(3), ...
    pars.probs(4), pars.boundary);

fprintf('%s\n', name);                
if ~ pars.update && exist([path, name, '.mat'], 'file')
    load([path, name, '.mat']);
    disp('file exists, loaded');
else
    disp('update or file not exist');
    
    %%
    simul = Simulator(varargin{:}, 'plotting', 0);
    simul.simulate();
    
    %% measure
    ls = 1: ceil(pars.cir / pars.num): pars.cir;
    entropies = zeros(size(ls));
    for idx = 1: length(ls)
        tableau = simul.tableau.clone();
        qubits2trace = simul.check_generator.len2qubits(ls(idx), 'start', 1);
        tableau.partial_trace(qubits2trace);
        entropies(idx) = tableau.get_entropy();
    end
    
    pars_tmp = pars;
    clear pars
    save([path, name, '.mat']);
    fprintf('data saved:\n %s\n', [path, name]);
    fprintf('Time elapsed: %f s\n',  etime(clock(), tictime));
    pars = pars_tmp;
end
    
output = {entropies};
xs = ls / pars.cir;
%%
if pars.plotting    
    title_str = {'title', 'sub'};
    path = sprintf('./figures_%s/', 'bipartite');
    mkdir(path);
    save_str = [path, name];  
    plot_bipartite(output, xs, title_str, save_str, ...
        'legend_str', {sprintf('L=%d', pars.cir)});
end

end
