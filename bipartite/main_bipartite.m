tictime = clock(); 
script_path = fileparts(mfilename('fullpath'));
addpath(genpath([script_path, '/../../utils']));
addpath([script_path, '/../']);

rng('shuffle');


%%
ps = flipud([0; 0.25; 0.683; 0.8]);

for i = 1:length(ps)
    p = ps(i);
    pxyz = (1 - p)/3;
    bipartite('probs', [pxyz, pxyz, pxyz, p], ...
        'cir', 12, 'update', 0);
end

bipartite('probs', [0.1, 0.1, 0.8, 0], ...
    'cir', 12, 'update', 0);
