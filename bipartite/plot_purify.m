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
