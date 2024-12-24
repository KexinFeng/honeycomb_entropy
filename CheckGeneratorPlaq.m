classdef CheckGeneratorPlaq < CheckGenerator
properties
    probs
end
methods   
    function obj = CheckGeneratorPlaq(varargin)
        obj@CheckGenerator(varargin{:});
        
        ip = inputParser;
        ip.KeepUnmatched = true;
        ip.addParameter('probs', [1/4, 1/4, 1/4, 1/4]);
        ip.parse(varargin{:});
        pars = ip.Results;

        obj.probs = pars.probs;
    end

    function set_probs(obj, probs)
        obj.probs = probs;
    end

    % function [row, check] = generate_bond(obj, x, y)
    %     cir = obj.cir;
    %     wid = obj.wid;
    % 
    %     % Generate bond
    %     if strcmp(obj.boundary, "periodic")
    %         check = randi([1, 4]); % 1:z, 2:x, 3:y, 4:plaq
    %     elseif strcmp(obj.boundary, 'open')
    %         if x==1 && y==1 
    %             check = randsample([1, 3], 1);
    %         elseif x == cir && y == wid
    %             check = [];
    %             row = [];
    %             return
    %         elseif x == 1
    %             check = randi([1, 3]);
    %         elseif y == wid
    %             check = randsample([2, 3], 1);
    %         elseif x == cir
    %             check = randsample([1, 2, 4], 1);
    %         else
    %             check = randi([1, 4]);
    %         end
    %     else
    %         error('%s not right', obj.boundary)
    %     end
    % 
    %     if check < 4
    %         row = obj.bond2row(x, y, check);
    %     else
    %         row = obj.plaq2row(x, y);
    %     end
    % 
    %     pauli = Util.row2pauli(row, cir, wid);
    %     fprintf('check %s %d %d: %s\n', obj.check_names{check}, x, y, pauli)
    % end


    function [row, check] = generate_bond(obj, x, y)
        probs = obj.probs; 

        % 1:z, 2:x, 3:y, 4:plaq
        cir = obj.cir;
        wid = obj.wid;
        
        assert(abs(sum(probs) - 1) < 1e-3);       
        
        function check = sample(probs)
            probs = probs / sum(probs);
            edges = [0, cumsum(probs)];
            check = find(edges < rand(), 1, 'last'); 
        end

        % Generate bond
        if strcmp(obj.boundary, "periodic")
            check = sample(probs);
        elseif strcmp(obj.boundary, 'open')
            if x==1 && y==1 
                % check = randsample([1, 3], 1);
                probs([2, 4]) = 0;
                check = sample(probs);
            elseif x == cir && y == wid
                check = [];
                row = [];
                return
            elseif x == 1
                % check = randi([1, 3]);
                probs(4) = 0;
                check = sample(probs); 
            elseif y == wid
                % check = randsample([2, 3], 1);
                probs([1, 4]) = 0;
                check = sample(probs);                
            elseif x == cir
                % check = randsample([1, 2, 4], 1);
                probs(3) = 0;
                check = sample(probs);
            else
                % check = randi([1, 4]);
                check = sample(probs);
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
        fprintf('check %s %d %d: %s\n', obj.check_names{check}, x, y, pauli)
    end


    %% Utility
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

end
end

