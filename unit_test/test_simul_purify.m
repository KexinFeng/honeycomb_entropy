filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
rng(24);

%% Generator (deprecated)
simul = Simulator('cir', 2, 'wid', 2, 'boundary', 'periodic');
generator = CheckGenerator('cir', 2, 'wid', 2, 'boundary', 'periodic');
simul.check_generator = generator;
simul.verbose = true;
simul.T = 20;
simul.simulate()

%% GeneratorPlaquet
simul = Simulator('cir', 2, 'wid', 2, 'boundary', 'periodic');
simul.check_generator.set_probs([1/3, 1/3, 1/3, 0])
simul.verbose = true;
simul.T = 20;
simul.simulate()

%% GeneratorPlaquet
rng(24);
simul = Simulator('cir', 2, 'wid', 2, 'boundary', 'periodic');
generator = CheckGeneratorPlaq('cir', 2, 'wid', 2, 'boundary', 'periodic');
simul.check_generator = generator;
simul.verbose = true;
simul.T = 20;
simul.simulate()

%%
generator = simul.check_generator;

generator.set_probs([1/3, 1/3, 0, 1/3]);
simul.simulate();

generator.set_probs([0, 0, 1/2, 1/2]);
simul.simulate();

generator.set_probs([1/3, 1/3, 0, 1/3]);
simul.simulate();

clc
fprintf('\n------------\n')
generator.set_probs([0, 0, 1/2, 1/2]); % z:1, x:2, y:3, p:4
simul.simulate();

dbstop = 1;

