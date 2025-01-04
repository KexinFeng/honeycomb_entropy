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
    
    % Convert pars to a cell array of parameter-value pairs
    fields = fieldnames(pars);      
    values = struct2cell(pars);     
    varargin = reshape([fields, values]', 1, []); 
    
    %% get late-time state
    res = purify(varargin{:});
    
    %% measure bipartite
    sys_sizes = 0: ceil(pars.cir / pars.num): pars.cir;
    sys_sizes = unique(sort([sys_sizes, floor(pars.cir / 2)]));
    
    entropies = zeros(size(sys_sizes));
    for idx = 1: length(sys_sizes)
        tableau = res.simul.tableau.clone();
        sys_size = sys_sizes(idx);
        qubits_env = res.simul.check_generator.len2qubits(pars.cir - sys_size, 'start', 1 + sys_size);
        tableau.partial_trace(qubits_env);
        entropies(idx) = tableau.get_entropy();
    end
    % Delta S
    idx = (sys_sizes == pars.cir / 2);
    delta_S = entropies - entropies(idx);
    
    output = {delta_S / pars.cir, entropies / pars.cir};
    xs = sys_sizes / pars.cir;
    
    
    %% plotting
    if pars.plotting
        name = sprintf('cir_%d_T_%d_probs_%.2f_%.2f_%.2f_%.2f_%s',...
        pars.cir, pars.T, pars.probs(1), pars.probs(2), pars.probs(3), ...
        pars.probs(4), pars.boundary);
    
        title_str = name;
        title_str = regexprep(title_str, '(?<=\D)_', '=');
        title_str = regexprep(title_str, '(?<=\d)_', ' ');
        
        fig_folder = [script_path, '/figures_bipartite/'];
        mkdir(fig_folder);
        save_str = [fig_folder, '/bipartite_', name];  
        plot_bipartite(output, xs, title_str, save_str, ...
            'legend_str', {sprintf('L=%d', pars.cir)});
    end

end


function plot_bipartite(outputs, xs, title_str, save_str, varargin)
    ip = inputParser;
    ip.addParameter('err', {});
    ip.addParameter('legend_str', {'L=10', 'pp', 'hh'});
    
    ip.parse(varargin{:})
    pars = ip.Results;
    
    %%
    ys = outputs{1};
    
    prep_plots
    figure
    hold on
    box on
    if isempty(pars.err)
        hb = plot(xs, ys, '*b-');
    else
        hb = errorbar(xs, ys, pars.err{1}, '*b-');
    end
    
    xlabel('l/L', 'Interpreter', 'latex')
    ylabel('S/Lln2', 'Interpreter', 'latex')
    
    legend(hb, pars.legend_str, 'AutoUpdate', 'off', 'Location', 'southeast');
    title(title_str);
    
    
    %%
    if ~isempty(save_str)
        str = save_str;
        print(gcf, '-depsc2', [str,'.eps']);
        eps2pdf([str,'.eps'], [str,'.pdf'], 1);
        delete([str,'.eps']);
    end
end


