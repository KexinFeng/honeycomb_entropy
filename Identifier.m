classdef Identifier < handle
properties
    cir
    wid
    shift
    % internal
    Ns
end

methods
    function obj = Identifier(varargin)
        ip = inputParser;
        ip.KeepUnmatched = true; 
        ip.PartialMatching = false; 
        ip.addParameter('cir', 4);
        ip.addParameter('wid', 4);
        ip.addParameter('shift', 0); % x shift per y increase. 

        % Parse input arguments
        ip.parse(varargin{:});
        pars = ip.Results;

        % Loop through the fields of pars and assign to properties
        fields = fieldnames(pars);
        for i = 1:numel(fields)
            field = fields{i};
            obj.(field) = pars.(field); 
        end
        
        % Internal
        obj.Ns = obj.cir * obj.wid * 2;
    end

    
    function [xs_out, ys_out] = comb(obj, xs_in, ys_in)
        if (isnumeric(obj.shift) && obj.shift == -1/2) || ...
                ((ischar(obj.shift) || isstring (obj.shift)) && strcmp(obj.shift, 'rect')) 
            if mod(obj.wid, 2) ~= 0
                error('Incorrect geometry')
            end
            [ys_out, quo] = homod(ys_in, obj.wid);  % homod(ys, wid) <=> ys -= wid
            x_shift = (quo - 1) * obj.wid / 2;
            xs_out = homod(xs_in + x_shift, obj.cir);

        elseif obj.shift == 0
            xs_out = homod(xs_in, obj.cir);
            ys_out = homod(ys_in, obj.wid);
            
        else
            if mod(obj.wid, 2) ~= 0
                error('Incorrect geometry')
            end
            [ys_out, quo] = homod(ys_in, obj.wid);
            x_shift = (quo - 1) * obj.wid * (-obj.shift);
            xs_out = homod(xs_in + x_shift, obj.cir);
        end
    end
end
end