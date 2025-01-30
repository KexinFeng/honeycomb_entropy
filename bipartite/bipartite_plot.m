function bipartite_plot(varargin)
tictime = clock(); 
script_path = fileparts(mfilename('fullpath'));
addpath(genpath([script_path, '/../../utils']));
addpath([script_path, '/../']);


%% input
ip = inputParser;
ip.KeepUnmatched = true;
ip.PartialMatching = false;

% ip.addParameter('cir', 14);
ip.addParameter('boundary', 'periodic');
ip.addParameter('shift', 0);
ip.addParameter('T', 150);
ip.addParameter('verbose', false);
% ip.addParameter('probs', [0.25, 0.25, 0.25, 0.25]);
ip.addParameter('num', 10);
ip.addParameter('plotting', 1);
ip.addParameter('update', 0);
ip.addParameter('init', 'flux_free');

ip.parse(varargin{:});
pars = ip.Results;

% Convert pars to a cell array of parameter-value pairs
fields = fieldnames(pars);      
values = struct2cell(pars);     
varargin = reshape([fields, values]', 1, []); 

%% Prepare input
cirs = [18, 24, 30, 36, 42, 48, 54, 60];  
cirs = [6, 10, 12, 14, 16, 18, 24];  
% cirs = [10];
probs = [0.25, 0.25, 0.25, 0.25];



%% save
name_file = sprintf('T_%d_probs_%.2f_%.2f_%.2f_%.2f_%s',...
    pars.T, probs(1), probs(2), probs(3), ...
    probs(4), pars.boundary);

%% Plotting
colors = linspecer(length(cirs), 'sequential');
% colors = cbrewer('qual', 'Paired', 6);
cmap = colormap(linspecer);

figure
prep_plots

hold on
legends = {};
hands = {};

for i = 1: length(cirs)
    L = cirs(i);

    % if L >= 36
    %     pars.T = 300;
    % end

    %% save
    name = sprintf('cir_%d_T_%d_probs_%.2f_%.2f_%.2f_%.2f_%s',...
        L, pars.T, probs(1), probs(2), probs(3), ...
        probs(4), pars.boundary);
    
    fprintf('%s\n', name);                
    data_folder = sprintf([script_path, '/data_%s/'], replace(mfilename, '_plot', ''));

    loaded = load([data_folder, name, '.mat']);


    h = plot(loaded.ls / L, loaded.entropies / L, '-', 'Color', colors(i, :));
    legends = [legends, sprintf('L=%d', L)];
    hands = [hands, h];

    lgd = legend(hands, legends, 'Location', 'best');
    set(lgd, 'AutoUpdate', 1);

    pause(0.001);
    
end        


title_str = name_file;
title_str = regexprep(title_str, '(?<=\D)_', '=');
title_str = regexprep(title_str, '(?<=\d)_', ' ');


%% Plot setting
grid off
box on

title(title_str)

lgd = legend(hands, legends, 'Location', 'south');
% legend('boxoff');
set(lgd, 'AutoUpdate', 1);

dbstop = 1;

%% Plot saving
plot_folder = [script_path, '/figures_bipartite_plot/'];
triexpf(~exist(plot_folder, 'dir'), {@mkdir, plot_folder}, {});
plot_file = sprintf('/probs_%.2f_%.2f_%.2f_%.2f_periodic',...
    probs(1), probs(2), probs(3), probs(4));

savepdf([plot_folder, plot_file]);
end

