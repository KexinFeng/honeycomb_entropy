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


legend(hb, pars.legend_str, 'AutoUpdate','off');
title(title_str);


%%
if ~isempty(save_str)
    str = save_str;
    print(gcf, '-depsc2', [str,'.eps']);
    eps2pdf([str,'.eps'], [str,'.pdf'], 1);
    delete([str,'.eps']);
end
end
