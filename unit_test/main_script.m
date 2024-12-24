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
simul.T = 10;
simul.simulate()

%%
generator = simul.check_generator;

generator.set_probs([1/3, 1/3, 0, 1/3]);
% simul.check_scenario
simul.simulate();

generator.set_probs([0, 0, 1/2, 1/2]);
simul.simulate();

%% GeneratorPlaquet
simul = Simulator('cir', 2, 'wid', 2, 'boundary', 'open');
simul.verbose = true;
simul.T = 10;
simul.simulate()

%%
generator = simul.check_generator;

generator.set_probs([1/3, 1/3, 0, 1/3]);
simul.simulate();

generator.set_probs([0, 0, 1/2, 1/2]);
simul.simulate();

dbstop = 1;

