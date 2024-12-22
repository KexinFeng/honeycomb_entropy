p = 2; % Characteristic of the field (prime number)
m = 1; % Degree of the field
X = [1, 0, 1; 
     1, 1, 1;
     1, 0, 0];

Xg = gf(X, m, p);
rank(Xg)

