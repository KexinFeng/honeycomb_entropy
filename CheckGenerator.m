classdef CheckGenerator < handle
properties
    cir
    wid
    boundary
    % internal
    Ns
    check_names
end

methods
    function obj = CheckGenerator(varargin)
        ip = inputParser;
        ip.KeepUnmatched = true; % Allow unmatched parameters
        ip.addParameter('cir', 2);
        ip.addParameter('wid', 2);
        ip.addParameter('boundary', 'open');

        % Parse input arguments
        ip.parse(varargin{:});
        pars = ip.Results;

        % Loop through the fields of pars and assign to properties
        fields = fieldnames(pars);
        for i = 1:numel(fields)
            field = fields{i};
            obj.(field) = pars.(field); % Dynamic field assignment
        end

        obj.Ns = obj.cir * obj.wid * 2;

        obj.check_names = {'Z', 'X', 'Y', 'P'};
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


    function [row, check] = generate_check(obj, x, y)
        cir = obj.cir;
        wid = obj.wid;

        % Generate bond
        if strcmp(obj.boundary, "periodic")
            check = randi([1, 4]); % 1:z, 2:x, 3:y, 4:plaq
        elseif strcmp(obj.boundary, 'open')
            if x==1 && y==1 
                check = randsample([1, 3], 1);
            elseif x == cir && y == wid
                check = [];
                row = [];
                return
            elseif y == wid
                check = randsample([2, 3], 1);
            elseif x == cir
                check = randsample([1, 2], 1);
            else
                check = randi([1, 4]);
            end
        else
            error('%s not right', obj.boundary)
        end

        if check < 4
            row = obj.bond2row(x, y, check);
        else
            row = obj.plaq2row(x, y);
        end

        pauli = Util.row2pauli(row, cir, wid);
        fprintf('%s\n check %d %d: %s\n', obj.check_names{check}, x, y, pauli)
    end


    %% Utilities
    function row = plaq2row(obj, x, y)
        cir = obj.cir;
        wid = obj.wid;
        Ns = obj.Ns;
    
        row = zeros(1, Ns*2);
        % (x, y, 0)-X-(x, y , 1)
        % (x-1, y+1, 0)-Z-(x-1, y , 1)
        % (x, y+1, 0)-Y-(x-1, y+1 , 1)
        
        % X
        bin = [1; 0];
        bin_arr = kron(bin, [1; 1]);
        xs = [x, x];
        ys = [y, y];
        ab = [0, 1]; 
        row = obj.fill_row_entry(row, bin_arr, xs, ys, ab);

        % Z
        bin = [0; 1];
        bin_arr = kron(bin, [1; 1]);
        xs = [homod(x-1, cir), homod(x-1, cir)];
        ys = [homod(y+1, wid), y];
        ab = [0, 1];
        row = obj.fill_row_entry(row, bin_arr, xs, ys, ab);

        % Y
        bin = [1; 1];
        bin_arr = kron(bin, [1; 1]);
        xs = [x, homod(x-1, cir)];
        ys = [homod(y+1, wid), homod(y+1, wid)];
        ab = [0, 1];
        row = obj.fill_row_entry(row, bin_arr, xs, ys, ab);

        % check_names = ['Z', 'X', 'Y'];
        % pauli = Util.row2pauli(row, cir, wid);
        % fprintf('%s\n pauli: %s\n', bond_names{bond}, pauli)
    end


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
        
        xs = [x, homod(x+(bond == 3), cir)];
        ys = [y, homod(y+(bond == 1), wid)];
        ab = [1, 0];
        
        row = obj.fill_row_entry(row, bin_arr, xs, ys, ab);
        % bond_names = ['Z', 'X', 'Y'];
        % pauli = Util.row2pauli(row, cir, wid);
        % fprintf('%s\n pauli: %s\n', bond_names{bond}, pauli)
    end


    function row = fill_row_entry(obj, row, bin_arr, xs, ys, ab)
        cir = obj.cir;
        wid = obj.wid;
        Ns = obj.Ns;
    
        % Site index: x = 1~cir; y = 1~wid; ab = 0,1
        pos = @(x, y, ab) y + wid.*(x-1) + cir*wid.*ab;
        idx = pos(xs, ys, ab);
        idx = [idx, idx + Ns];
        row(idx) = reshape(bin_arr, size(idx)); 
    end

end
end

