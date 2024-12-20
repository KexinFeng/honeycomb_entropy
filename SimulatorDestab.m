classdef SimulatorDestab
methods
    function obj = SimulatorDestab(inputArg1,inputArg2)
        %SIMULATORDESTAB Construct an instance of this class
        %   Detailed explanation goes here
        obj.Property1 = inputArg1 + inputArg2;
    end
    
    function outputArg = measure_simul(obj, inputArg)
        outputArg = obj.Property1 + inputArg;
    end
end
end

