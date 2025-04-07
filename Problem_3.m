% Double integral of f(x, y) = 2y*sin(x) + cos^2(x)
% over y = sin(x) to y = cos(x), x = 0 to pi/4

clc; clear; close all;

f = @(x, y) 2 .* y .* sin(x) + cos(x).^2;

% Exact value (symbolic integration)
syms x y
f_sym = 2*y*sin(x) + cos(x)^2;
I_exact = double(int(int(f_sym, y, sin(x), cos(x)), x, 0, pi/4));

% (a) Simpson's Rule with n = m = 4
n = 4; m = 4;
a = 0; b = pi/4;
hx = (b - a)/(2*n);
x = linspace(a, b, 2*n+1);
Isimp = 0;

for i = 1:2*n+1
    xi = x(i);
    yi_low = sin(xi);
    yi_high = cos(xi);
    hy = (yi_high - yi_low)/(2*m);
    y = linspace(yi_low, yi_high, 2*m+1);

    for j = 1:2*m+1
        weight_y = 1;
        if mod(j-1, 2) == 1
            weight_y = 4;
        elseif j ~= 1 && j ~= 2*m+1
            weight_y = 2;
        end
        fy(j) = weight_y * f(xi, y(j));
    end
    inner = hy/3 * sum(fy);

    weight_x = 1;
    if mod(i-1, 2) == 1
        weight_x = 4;
    elseif i ~= 1 && i ~= 2*n+1
        weight_x = 2;
    end
    Isimp = Isimp + weight_x * inner;
end

Isimp = hx/3 * Isimp;

% (b) Gaussian Quadrature with n = m = 3
xi = [-sqrt(3/5), 0, sqrt(3/5)];
wi = [5/9, 8/9, 5/9];
Igauss = 0;

for i = 1:3
    x_i = 0.5*(b - a)*xi(i) + 0.5*(b + a);
    Jx = (b - a)/2;
    ylow = sin(x_i);
    yhigh = cos(x_i);
    Jy = (yhigh - ylow)/2;
    sum_y = 0;
    for j = 1:3
        y_j = 0.5*(yhigh - ylow)*xi(j) + 0.5*(yhigh + ylow);
        sum_y = sum_y + wi(j) * f(x_i, y_j);
    end
    Igauss = Igauss + wi(i) * Jy * sum_y * Jx;
end

% Output Results
fprintf('Exact Value            = %.6f\n', I_exact);
fprintf("a. Simpson's Rule      = %.6f\n", Isimp);
fprintf('b. Gaussian Quadrature = %.6f\n', Igauss);

% Compare these results with the exact value
eSimp_abs = abs(Isimp - I_exact);
eSimp_rel = abs(eSimp_abs/I_exact);
eGauss_abs = abs(Igauss - I_exact);
eGauss_rel = abs(eGauss_abs/I_exact);
fprintf("c. Simpson's Rule Absolute Error      = %.6f, Relative Error = %.4f%%\n", eSimp_abs, eSimp_rel*100);
fprintf("   Gaussian Quadrature Absolute Error = %.6f, Relative Error = %.4f%%\n", eGauss_abs, eGauss_rel*100);