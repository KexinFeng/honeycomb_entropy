function [row, bond, res] = generate_bond(x, y, pars)
    cir = pars.cir;
    wid = pars.wid;
    
    % Generate bond
    if strcmp(pars.boundary, "periodic")
        bond = randi([1, 3]); % 1:z, 2:x, 3:y
    elseif strcmp(pars.boundary, 'open')
        if x==1 && y==1 
            bond = randsample([1, 3], 1);
        elseif x == cir && y == wid
            bond = [];
            return
        elseif y == wid
            bond = randsample([2, 3], 1);
        elseif x == cir
            bond = randsample([1, 2], 1);
        else
            bond = randi([1, 3]);
        end
    else
        error('%s not right', boundary)
    end
    
    row = bond2row(x, y, bond, pars);

    res.bond2row = @bond2row;
end

function row = bond2row(x, y, bond, pars)
    % bond_names = ['Z', 'X', 'Y'];

    cir = pars.cir;
    wid = pars.wid;
    Ns = cir * wid * 2;

    row = zeros(1, Ns*2);
    % Site index: x = 1~cir; y = 1~wid; ab = 0,1
    pos = @(x, y, ab) y + wid.*(x-1) + cir*wid.*ab;

    % (x, y, 0)-X-(x, y , 1)
    % (x, y+1, 0)-Z-(x, y, 1)
    % (x+1, y, 0)-Y-(x, y, 1)
    bin = dec2bin(bond, 2);
    bin_arr = kron(str2num(bin(:)), [1; 1]);
    
    xs = [x, homod(x+(bond == 3), cir)];
    ys = [y, homod(y+(bond == 1), wid)];
    ab = [1, 0];
    
    idx = pos(xs, ys, ab);
    idx = [idx, idx + Ns];
    
    row(idx) = reshape(bin_arr, size(idx));
    % bond_names = ['Z', 'X', 'Y'];
    % pauli = Util.row2pauli(row, cir, wid);
    % fprintf('%s\n pauli: %s\n', bond_names{bond}, pauli)
end

