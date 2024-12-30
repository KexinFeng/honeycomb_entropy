filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
rng(24);


%% GeneratorPlaquet
simul = Simulator('cir', 10, 'boundary', 'periodic');
simul.check_generator.set_probs([0, 1/5, 1/3, 2/3-1/5]);
simul.check_generator.set_probs([1/4, 1/4, 1/4, 1/4]);
simul.verbose = false;
simul.T = 200;
simul.simulate()


dbstop = 1;

