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
    cirs = [6, 10, 12, 14, 16, 18, 24, 30, 36, 42];  
    % cirs = [10];
    probs = [0.25, 0.25, 0.25, 0.25];

    %% Pre-plot purify
    for i = 1: length(cirs)
        L = cirs(i);
    
        ls = 0: max(1, floor(L / pars.num)): L;
        ls = unique(sort([ls, floor(L / 2)]));
            
        % Compute entropy
        res = purify(varargin{:}, 'cir', L, 'probs', probs, 'plotting', 1);

        pause(0.001);
    end
    
    
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
    
        ls = 0: max(1, floor(L / pars.num)): L;
        ls = unique(sort([ls, floor(L / 2)]));
            
        % Compute entropy
        res = purify(varargin{:}, 'cir', L, 'probs', probs, 'plotting', 0);
        entropies = zeros(size(ls));

        for idx = 1: length(ls)
            tableau = res.simul.tableau.clone();
            sys_size = ls(idx);
            qubits_env = res.simul.check_generator.len2qubits(L - sys_size, 'start', 1 + sys_size);
            tableau.partial_trace(qubits_env);
            entropies(idx) = tableau.get_entropy();
        end  

        % Delta S
        % idx = (ls == L / 2);
        % delta_S = entropies - entropies(idx);
    
        h = plot(ls / L, entropies / L, '-', 'Color', colors(i, :));
        legends = [legends, sprintf('L=%d', L)];
        hands = [hands, h];

        lgd = legend(hands, legends, 'Location', 'best');
        set(lgd, 'AutoUpdate', 1);

        pause(0.001);
    end
        
    
    %% Plot setting
    grid off
    box on
    
    lgd = legend(hands, legends, 'Location', 'south');
    % legend('boxoff');
    set(lgd, 'AutoUpdate', 1);
    
    dbstop = 1;

    %% Plot saving
    save_str = [script_path, ...
        sprintf('/figures_bipartite_plot/probs_%.2f_%.2f_%.2f_%.2f',...
        probs(1), probs(2), probs(3), probs(4))];
    if ~isempty(save_str)
        str = save_str;
        print(gcf, '-depsc2', [str,'.eps']);
        eps2pdf([str,'.eps'], [str,'.pdf'], 1);
        delete([str,'.eps']);
    end
end

% 
% 
% 
%     % Add legend
%     legend_entries = arrayfun(@(L) sprintf('L = %d', L), L_values, 'UniformOutput', false);
%     legend(legend_entries, 'Location', 'best');
%     hold off;
% 
% 
% 
%     %% get late-time state
%     [output, xs] = bipartite(varargin{:}, 'cir', L);
% 
% 
% 
% 
%     %% data loading
%     tic
%     Ts_cand = linspace(-2, 1.5, 15);
%     Ts_cand = 10.^Ts_cand;
% 
%     Ts = [0.01, 0.05, 0.1, 0.5, 1,2, 3, 5, 10, 50];
%     Ts = [0.03, 0.05, 0.1, 0.5, 1,5, 10, 50];  
% 
% 
% %     if kap > 0
%     T = Ts(ord);
%     idx = find(Ts_cand >= T, 1);
%     T = Ts_cand(idx);
% %     else
% %         T = Ts(ord);
% %     end
% 
%     fprintf('T = %f\n', T);    
% 
%     colors = linspecer(10, 'sequential');
%     colors = cbrewer('qual', 'Paired', 6);
%     cmap = colormap(linspecer);
% 
%     prep_plots
% %     figure('Position', [91.0000  215.6667  486.6667  370.0000]);
% %     set(0,'defaultLineLineWidth',2.5);   % set the default line width to lw
% 
%     hold on
%     legends = {};
%     hands = {};
% 
%     output = straSamp_integ_kap_load_2('fun', fun, 'T', T, 'cir', cir, 'kap', kap, 'update', 0,...
%         'use_gpu', use_gpu, 'Nsample', Nsample, 'simDim', simDim);   
%     Ihats = output{1};
%     Ierrs = output{2};
%     Omgs = output{3};
%     nbar = output{4};
% 
% 
%     legend_strs = {'11', '12', '21', '22'};
%     for idx = 1:4
%         i = idx + triexp(strcmp(IR, 'real'), 4, 0);
%         c_idx = homod(idx + 1, 4);
% 
%         simDim =  length(Ihats{i});
%         step = 3;
%         idxs = [1:8, 9:step:simDim];
% 
%         x = Omgs(idxs);
%         y =  reshape(Ihats{i}, 1, []);
%         if strcmp(IR, 'real')
%             y = -y;
%         end
%         y = y(idxs);
%         yerr = reshape(Ierrs{i}, 1, []);
%         yerr = yerr(idxs);
% 
% %         plot(x, y, 'o');
% 
%         %     pp = spline(x, y);
%         method ='pchip';
%         xx = linspace(min(x), max(x), 200);
%         yy = interp1(x, y, xx, method);
%         h = plot(xx, yy, 'Color', colors(c_idx, :));
% 
% %         errorbar(x(idxs), y(idxs), yerr(idxs)./2, '*', 'Color', colors(c_idx, :), ...
% %             'MarkerSize', msz, 'LineWidth', lw-.5);
% 
%         legends = [legends, sprintf('%s', legend_strs{idx})];
%         hands = [hands, h];
%     end
% 
% 
% 
% %% plot setting
% %     legends = fliplr(legends);
% %     hands = fliplr(hands);
%     if strcmp(IR, 'imag')
%         ymax = 5;
%         axis([0, 12, 0, 4]);
%     else
%         ymax = 5;
%         axis([0, 12, -4, 2]);
%     end
%             xlim([0, 12]);
% 
% 
%     set(gca,'XTick',[0 3 6 9 12]);
%     if find(ord == [4, 5, 6])
%         xlabel('\Omega');
%     end
%     if find(ord == [1, 4])
%         if strcmp(IR, 'imag')
%             ylabel('$-$Im$\Pi_{mm^\prime}$', 'interpreter', 'latex');
%         else
%             ylabel('Re$\Pi_{mm^\prime}$', 'interpreter', 'latex');
%         end
%     end
%     if find(ord == [1,2,3])
%         set(gca,'xticklabel',[])
%     end
%     if find(ord == [2, 3, 5, 6])
%         set(gca,'yticklabel',[])
%     end
%     str_ord = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'};
%     text(gca, 0.02, 1.07 , sprintf('(%c)', str_ord{ord}),'Unit', 'normalized',  'FontSize', 14);
% 
% 
%     grid off
%     box on
% 
%     lgd = legend(hands, legends, 'Location', triexp(strcmp(IR, 'imag'), 'northeast', 'southeast'));
%     legend('boxoff');
%     set(lgd, 'AutoUpdate', 0);
% %     legend('NumColumns', 2);   
% 
% % %%
% %% legend
% string = {sprintf('T=%.2f, n=%.2f', T, nbar)};
% 
% ymax = ylim;
% ymax = ymax(2);
% 
% ys = [ymax*0.9, ymax*0.98];
% xs = [0.3, 0.4];
% 
% [xs, ys] = normal_pos(xs, ys);
% annotation(gcf,'textbox',...
%     [xs(1), ys(1), xs(2)-xs(1), ys(2)-ys(1)],...
%     'String',string,...
%     'Interpreter','tex',...
%     'FitBoxToText','on',...
%     'EdgeColor',[0, 0,0],...
%     'BackgroundColor',[1 1 1],...
%     'FontSize', 12,...
%     'LineWidth', 0.5);
% % 
% %%
% eps1 = 7.38-0.03;
% eps2 = 10.7-0.59;
% xstar= eps1;
% line([xstar, xstar], ylim,'LineStyle', ':', 'LineWidth', 1, 'Color', colors(end-1, :));
% xstar= eps2;
% line([xstar, xstar], ylim,'LineStyle', ':', 'LineWidth', 1, 'Color', colors(end, :));
% 
% %% gap
% if kap > 0
% wid = cir;
% H = hamiltonianT(cir, wid, 0, 'change_gauge', -1,'kap', kap);
% s = diagnlz2(H, [cir, wid], kap, 'use_gpu', 0);
% s = double(gather(s));
% xstar = 2*min(6*sqrt(3)*kap, min(2*(s))); % \Delta_kap
% % line([xstar, xstar], ylim, 'LineWidth', 0.5);
% plot(xstar, 0, '^', 'Color', colors(end, :));
% end
% %%
% 
% 
% %     fun_info = functions(fun);
% %     fun_name = fun_info.function;
% %     fig_path = sprintf('./figures_%s/', fun_name);
% %     mkdir(fig_path);
% %     save_str = sprintf('IntenMc_%d_T_%.2f', cir, T);
% %     savepdf([fig_path, save_str])
% 
% end