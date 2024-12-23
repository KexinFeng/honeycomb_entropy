filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
rng(24);

%%
simul = Simulator('cir', 2, 'wid', 2, 'boundary', 'open');
simul.verbose = true;
simul.T = 10;
simul.simulate()


%% GeneratorPlaquet
simul = Simulator('cir', 2, 'wid', 2, 'boundary', 'open');
generator = CheckGeneratorPlaq('cir', 2, 'wid', 2, 'boundary', 'open');
simul.check_generator = generator;
simul.verbose = true;
simul.T = 20;
simul.simulate()
