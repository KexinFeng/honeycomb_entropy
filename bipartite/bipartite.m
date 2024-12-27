function [output, xs] = bipartite(varargin)
tic 
filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../../utils']));
addpath([filepath, '/../']);

rng(24);
clc;

%% input
ip = inputParser;
ip.addParameter('cir', 10);
ip.addParameter('wid', 10);
ip.addParameter('plotting', 1);
ip.addParameter('boundary', 'periodic');
ip.addParameter('T', 30);
ip.addParameter('verbose', false);
ip.addParameter('probs', [1/4, 1/4, 1/4, 1/4]);
ip.addParameter('num', 8);

ip.parse(varargin{:});
pars = ip.Results;

pars.wid = cir;


%%
simul = Simulator(varargin{:});
simul.simulate();

%% measure
ls = 1: ceil(pars.cir / num): pars.cir;
entropies = zerosLike(ls);
for idx = 1:length(ls)
    tableau = simul.tableau.clone();
    qubits = []; % f(l)
    tableau.partial_trace(qubits);
    entropies(idx) = tableau.get_entropy();
end


output = {entropies};
xs = ls;

end
