simul = Simulator('cir', 2, 'wid', 2, 'boundary', 'open');
simul.verbose = true;
simul.T = 30;
simul.simulate()