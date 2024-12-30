filepath = fileparts(mfilename('fullpath'));
addpath(genpath([filepath, '/../']));
addpath(genpath([filepath, '/../../utils']));
rng(24);

%%
simul = Simulator('cir', 2, 'wid', 2, 'boundary', 'periodic', 'verbose', 1);
cir = simul.cir;
wid = simul.wid;
Ns = cir * wid *2;

tableau = simul.tableau;
check_generator = simul.check_generator;

%% init
simul.zero_flux_init()


%% 
dbstop = 1;
disp('success')
