function res = bipartite(varargin)
tic 
filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../utils']));
% addpath([filepath, '/util.m']);

rng(24);
clc;
res = struct();

%% input
ip = inputParser;
ip.addParameter('cir', 10);
ip.addParameter('wid', 10);
ip.addParameter('plotting', 1);
ip.addParameter('boundary', 'open');
ip.addParameter('T', 30);
ip.addParameter('testing', false);
ip.addParameter('verbose', false);

ip.parse(varargin{:});
pars = ip.Results;








end
