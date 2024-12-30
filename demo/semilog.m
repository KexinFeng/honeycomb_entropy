x = 1:10;
y = 10.^x;

figure
% Logarithmic plot with log scale on y-axis
semilogy(x, y);
grid on;
xlabel('x');
ylabel('log(y)');
title('Logarithmic Scale on Y-Axis');
