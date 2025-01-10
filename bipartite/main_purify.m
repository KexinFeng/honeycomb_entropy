function main_purify(varargin)
    tictime = clock(); 
    script_path = fileparts(mfilename('fullpath'));
    addpath(genpath([script_path, '/../../utils']));
    addpath([script_path, '/../']);
    
    rng('shuffle');
    %% input
    ip = inputParser;
    ip.KeepUnmatched = true;
    ip.PartialMatching = false;
    
    % ip.addParameter('cir', 10);
    ip.addParameter('boundary', 'periodic');
    ip.addParameter('shift', 0);
    ip.addParameter('T', 300);
    ip.addParameter('verbose', false);
    % ip.addParameter('probs', [0.25, 0.25, 0.25, 0.25]);
    ip.addParameter('num', 10);
    ip.addParameter('plotting', 0);
    ip.addParameter('update', 0);
    ip.addParameter('init', 'flux_free');
    % ip.addParameter('use_gpu', false);
    
    
    ip.parse(varargin{:});
    pars = ip.Results;
    
    % Convert pars to a cell array of parameter-value pairs
    fields = fieldnames(pars);      
    values = struct2cell(pars);     
    varargin = reshape([fields, values]', 1, []); 
    

    %% prepare inputs
    ps = flipud([0; 0.25; 0.683; 0.8]);
    probs = zeros(length(ps), 4);
    for i = 1: length(ps)
        p = ps(i);
        pxyz = (1 - p)/3;
        probs(i, :) = [pxyz, pxyz, pxyz, p];
    end
    % probs = [probs; [0.1, 0.1, 0.8, 0]];
    % probs_idx = 1: size(probs, 1);
    % cirs = [18, 24, 30, 36, 42, 48, 54, 60];
    
    % Trial
    cirs = [10, 12, 18, 24, 30, 36, 42, 48, 54, 60];  
    cirs = [14, 24, 42, 48, 54, 60];
    cirs = [18, 24, 30, 36, 42];
    
    cirs = [18, 24, 48, 54, 60];
    cirs = [36, 42, 48, 54, 60];
    % pars.T = 20;

    probs = [0.25, 0.25, 0.25, 0.25];
    probs_idx = 1: size(probs, 1);

    % Cartesian product
    [cirs_x, probs_idx_y] = meshgrid(cirs, probs_idx);  % y-first-order
    overheads = zeros(size(cirs_x));

    %% 
    %launch_parpool('spare_numC', 0);
    %parfor ord = 1: length(cirs_x)
    for ord = 1: length(cirs_x)

        tictime1 = clock();
        purify(varargin{:}, ...
            'probs', probs(probs_idx_y(ord), :), ...
            'cir', cirs_x(ord), ...
            'update', pars.update, ...
            'plotting', 0, ...
            'T', pars.T);
        
        overheads(ord) = etime(clock(), tictime1);
        fprintf('cir=%d, Time elapsed: %f s\n', cirs_x(ord), etime(clock(), tictime1));
    end
    
    fprintf('Time elapsed: %f s\n', etime(clock(), tictime));

    overhead = mean(overheads, 1);
    
    disp(cirs_x(1, :))
    disp(overhead)

    %% save
    data_folder = sprintf([script_path, '/data_%s/'], mfilename);
    triexpf(~exist(data_folder, 'dir'), {@mkdir, data_folder}, {});
    name = sprintf('T_%d_cir-min_%d_cir-max_%d', pars.T, min(cirs), max(cirs));

    save([data_folder, name, '.mat']);

    %% plotting
    figure
    plot(cirs_x(1, :), overhead)
    xlim([0, 60])

    dbstop = 1;
    
    % %% plotting
    % % title_str = {'title', 'sub'};    
    % title_str = name;
    % title_str = regexprep(title_str, '(?<=\D)_', '=');
    % title_str = regexprep(title_str, '(?<=\d)_', ' ');
    % 
    % fig_folder = [script_path, '/figures_bipartite/'];
    % mkdir(fig_folder);
    % save_str = [fig_folder, '/bipartite_', name];  
    % plot_bipartite(output, xs, title_str, save_str, ...
    %     'legend_str', {sprintf('L=%d', pars.cir)});

end




