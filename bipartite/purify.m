function res = purify(varargin)
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
ip.addParameter('init', 'flux_free');

ip.parse(varargin{:});
pars = ip.Results;

pars.wid = pars.cir;
% Convert pars to a cell array of parameter-value pairs
fields = fieldnames(pars);      
values = struct2cell(pars);     
varargin = reshape([fields, values]', 1, []); 

%%
data_folder = sprintf([script_path, '/data_%s/'], mfilename);
triexpf(~exist(data_folder, 'dir'), {@mkdir, data_folder}, {});

name = sprintf('cir_%d_T_%d_probs_%.2f_%.2f_%.2f_%.2f_%s',...
    pars.cir, pars.T, pars.probs(1), pars.probs(2), pars.probs(3), ...
    pars.probs(4), pars.boundary);

fprintf('%s\n', name);                
if ~ pars.update && exist([data_folder, name, '.mat'], 'file')
    load([data_folder, name, '.mat']);
    disp('file exists, loaded');
else
    disp('update or file not exist');
    
    %%
    simul = Simulator(varargin{:});
    
    if strcmp(pars.init, 'flux_free')
        % init to be a single flux sector
        [es0, ts0] = simul.zero_flux_init();
        % simulate
        [es, ts] = simul.simulate();
        
        es = [es0, es];
        ts = [ts0, ts];
    else
        % simulate
        [es, ts] = simul.simulate();
        es = [simul.tableau.get_entropy(), es];
        ts = [0, ts];
    end

    %% save
    pars_tmp = pars;
    clear pars
    save([data_folder, name, '.mat']);
    fprintf('data saved:\n %s\n', [data_folder, name]);
    fprintf('Time elapsed: %f s\n',  etime(clock(), tictime));
    pars = pars_tmp;
end


%% plotting
if pars.plotting
    % title_str = {'title', 'sub'};    
    title_str = name;
    title_str = regexprep(title_str, '(?<=\D)_', '=');
    title_str = regexprep(title_str, '(?<=\d)_', ' ');
    
    fig_folder = [script_path, sprintf('/figures_%s/', mfilename)];
    mkdir(fig_folder);
    save_str = [fig_folder, 'purify_', name];  
    plot_purify(es, ts, title_str, save_str, ...
        'legend_str', {sprintf('L=%d', pars.cir)});
end

res.es = es;
res.ts = ts;
res.simul = simul;

end


function plot_purify(ys, xs, title_str, save_str, varargin)
ip = inputParser;
ip.addParameter('err', {});
ip.addParameter('legend_str', {'L=10', 'pp', 'hh'});

ip.parse(varargin{:})
pars = ip.Results;

%%
figure
prep_plots

hold on
box on

yaxis_min = log10(0.5);

hb = plot(xs, max(log10(ys), yaxis_min), '*b-');
grid on;

% Define the custom y-tick positions and labels
ytick_positions = log10([1, 2, 4, 10, 100]); 
ytick_labels = {'1', '2', '4', '10', '100'}; 
yticks(ytick_positions);
yticklabels(ytick_labels);

% Label
xlabel('time', 'Interpreter', 'latex')
ylabel('S/ln2', 'Interpreter', 'latex')

ylim([yaxis_min, log10(100)]);
xlim([0, 150]);
% axis([0, 125, 0, inf]);

legend(hb, pars.legend_str, 'AutoUpdate', 'off', 'Location', 'southeast');
title(title_str);


%% save
if ~isempty(save_str)
    str = save_str;
    print(gcf, '-depsc2', [str,'.eps']);
    eps2pdf([str,'.eps'], [str,'.pdf'], 1);
    delete([str,'.eps']);
end
end
