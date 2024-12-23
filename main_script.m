simul = Simulator('cir', 2, 'wid', 2);
simul.verbose = true;
simul.boundary = 'periodic';
simul.T = 10;
simul.simulate()