classdef CheckGenerator < handle
properties
    cir
    wid
    boundary
    % internal
    Ns
    check_names
    identifier
end

methods
    function obj = CheckGenerator(varargin)
        ip = inputParser;
        ip.KeepUnmatched = true; 
        ip.PartialMatching = false;
        ip.addParameter('cir', 4);
        ip.addParameter('wid', 4);
        ip.addParameter('boundary', 'periodic');

        % Parse input arguments
        ip.parse(varargin{:});
        pars = ip.Results;

        % Loop through the fields of pars and assign to properties
        fields = fieldnames(pars);
        for i = 1:numel(fields)
            field = fields{i};
            obj.(field) = pars.(field); % Dynamic field assignment
        end
        
        % Internal
        obj.Ns = obj.cir * obj.wid * 2;
        obj.check_names = {'Z', 'X', 'Y', 'P'};
        obj.identifier = Identifier(varargin{:});
    end

        
    %% Random generation
    function [row, bond] = generate_bond(obj, x, y)
        cir = obj.cir;
        wid = obj.wid;
        
        % Generate bond
        if strcmp(obj.boundary, "periodic")
            bond = randi([1, 3]); % 1:z, 2:x, 3:y
        elseif strcmp(obj.boundary, 'open')
            if x==1 && y==1 
                bond = randsample([1, 3], 1);
            elseif x == cir && y == wid
                bond = [];
                row = [];
                return
            elseif y == wid
                bond = randsample([2, 3], 1);
            elseif x == cir
                bond = randsample([1, 2], 1);
            else
                bond = randi([1, 3]);
            end
        else
            error('%s not right', obj.boundary)
        end
        
        row = obj.bond2row(x, y, bond);
    end


    %% Utilities
    function row = bond2row(obj, x, y, bond)
        % bond_names = ['Z', 'X', 'Y'];
        cir = obj.cir;
        wid = obj.wid;
        Ns = obj.Ns;
    
        row = zeros(1, Ns*2);

        % (x, y, 0)-X-(x, y , 1)
        % (x, y+1, 0)-Z-(x, y, 1)
        % (x+1, y, 0)-Y-(x, y, 1)
        bin = dec2bin(bond, 2);
        bin_arr = kron(str2num(bin(:)), [1; 1]);
        
        xs = [x+(bond == 3), x];
        ys = [y+(bond == 1), y];
        ab = [0, 1];
        [xs, ys] = obj.identifier.comb(xs, ys);
        row = obj.fill_row_entry(row, bin_arr, xs, ys, ab);
        % bond_names = ['Z', 'X', 'Y'];
        % pauli = Util.row2pauli(row, cir, wid);
        % fprintf('%s\n pauli: %s\n', bond_names{bond}, pauli)
    end
    

    function ret = pos(obj, x, y, ab)
        % Site index: x = 1~cir; y = 1~wid; ab = 0,1
        % pos = @(x, y, ab) y + wid.*(x-1) + cir*wid.*ab;
        ret = y + obj.wid .* (x-1) + obj.cir * obj.wid .* ab;
    end


    function [xs, ys, as] = pos_reverse(obj, lindex)
        % pos = @(x, y, ab) y + wid.*(x-1) + cir*wid.*ab;
        [ys, xs, as] = ind2sub([obj.wid, obj.cir, 2], lindex);
        % [res, a1] = homod(ind, cir*wid);
        % [res, x1] = homod(res, wid);
        % y1 = res;
        % assert(x==x1 && y==y1 && a==a1);
    end


    function row = fill_row_entry(obj, row, bin_arr, xs, ys, ab)
        % bin_arr represents a bilinear pauli operator @ (xs, ys, ab)
        idx = obj.pos(xs, ys, ab);
        idx = [idx, idx + obj.Ns];
        row(idx) = reshape(bin_arr, size(idx)); 
    end


    function ret = len2qubits(obj, delta_x, varargin)
        % [0, delta_x) right end excluded
        ip = inputParser;
        ip.KeepUnmatched = true; 
        ip.PartialMatching = false; 
        ip.addParameter('start', 1);

        % Parse input arguments
        ip.parse(varargin{:});
        pars = ip.Results;

        if ~strcmp(obj.boundary, 'periodic')
            error('Unsupported')
        end
        
        % delta's here mean the displacement from a reference point
        % At certain position along the circumference coordinate, a unit
        % cell contains two columns of sites, denoted by 1 and 2
        % respectively. Each column contains both A and B site.
        delta_ys_a_1 = 0: 2: obj.wid - 2;
        delta_xs_a_1 = 0: -1: -(obj.wid/2 - 1);
        
        delta_ys_b_1 = 1: 2: obj.wid - 1;
        delta_xs_b_1 = -1: -1: -obj.wid/2;
        
        delta_ys_a_2 = 1: 2: obj.wid - 1;
        delta_xs_a_2 = 0: -1: -(obj.wid/2 - 1);

        delta_ys_b_2 = 0: 2: obj.wid - 2;
        delta_xs_b_2 = 0: -1: -(obj.wid/2 - 1);

        unit_ys = [delta_ys_a_1, delta_ys_b_1, delta_ys_a_2, delta_ys_b_2];
        unit_xs = [delta_xs_a_1, delta_xs_b_1, delta_xs_a_2, delta_xs_b_2];
        unit_as = [zeros(size(delta_ys_a_1)), ones(size(delta_ys_b_1)), ...
            zeros(size(delta_ys_a_1)), ones(size(delta_ys_b_1))];

        translate_x = 0: delta_x - 1;
        xs = pars.start + unit_xs' + translate_x;
        ys = 1 + unit_ys' + zeros(size(translate_x));
        as = unit_as' + zeros(size(translate_x));
        
        [xs_rect, ys_rect] = obj.identifier.comb(reshape(xs, 1, []), reshape(ys, 1, []));
        as_rect = reshape(as, 1, []);
        ret = obj.pos(xs_rect, ys_rect, as_rect);
    end
end
end

