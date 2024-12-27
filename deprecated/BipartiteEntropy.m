classdef BipartiteEntropy  
properties
    cir
    wid
    boundary
    plotting
    T
    % internal
    simulator
end

methods
    function obj = BipartiteEntropy(varargin)
        filepath = fileparts(mfilename('fullpath'));
        addpath([filepath, '/util.m']);
        
        rng(24);
        clc;
    
        %% period boundary
        ip = inputParser;
        ip.KeepUnmatched = true;
        ip.addParameter('cir', 10);
        ip.addParameter('wid', 10);
        ip.addParameter('plotting', 1);
        ip.addParameter('boundary', 'open');
        ip.addParameter('T', 30);
        ip.addParameter('testing', false);
        ip.addParameter('verbose', false);
        
        ip.parse(varargin{:});
        pars = ip.Results;
        
        % Loop through the fields of pars and assign to properties
        fields = fieldnames(pars);
        for i = 1:numel(fields)
            field = fields{i};
            obj.(field) = pars.(field);  % Dynamic field assignment
        end

        obj.simulator = Simulator(varargin{:});
    end

        
    function res = runner(~, L)
        res = struct();
        kv_list = get_param(obj, 'cir', L, 'wid', L);
        simul = Simulator(kv_list{:});
        ls = 1: ceil(L/num): L;

        simul.simulate();
        
        % measure
        entropies = zerosLike(ls);
        for idx = 1:length(ls)
            tableau = simul.tableau.clone();
            qubits = []; % f(l)
            tableau.partial_trace(qubits);
            entropies(idx) = tableau.get_entropy();
        end
        
        res.entropies = entropies;
    end
end
end

