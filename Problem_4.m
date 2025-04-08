clc; clear; close all;

% ========== a) ∫ x^(-1/4) * sin(x) dx from 0 to 1 ==========
n = 4;
a = 0;
b = 1;
h = (b - a) / (2 * n);

x_vals = a : h : b;
G_vals = zeros(1, length(x_vals));

for i = 1:length(x_vals)
    G_vals(i) = G_func(x_vals(i));
end

% Composite Simpson’s Rule
I_G = h/3 * (G_vals(1) + 2*sum(G_vals(3:2:end-2)) + 4*sum(G_vals(2:2:end)) + G_vals(end));

% Analytic part of polynomial
syms x;
P_a = x^(-1/4) * (x - x^3/factorial(3));
I_P_a = double(int(P_a, 0, 1));

I_total_a = I_G + I_P_a;
fprintf('(a) Approximate value of ∫ x^{-1/4} * sin(x) dx from 0 to 1: %.7f\n', I_total_a);


% ========== b) ∫ x^(-4) * sin(x) dx from 1 to ∞ ==========
n = 4;              % Composite Simpson's rule n=4
a = 0;
b = 1;
h = (b - a) / (2*n);

t = linspace(a, b, 2*n+1);
f = zeros(1, 2*n+1);

for i = 1:length(t)
    if t(i) == 0
        f(i) = 0;
    else
        f(i) = t(i)^2 * sin(1 / t(i));
    end
end

% Composite Simpson’s Rule
I = (h / 3) * (f(1) + 4 * sum(f(2:2:end-1)) + 2 * sum(f(3:2:end-2)) + f(end));

fprintf('(b) Approximate value of ∫ x^{-4} * sin(x) dx from 1 to ∞: %.7f\n', I);
% ========== Function Definitions ==========

function y = G_func(x)
    if x == 0
        y = 0;
    elseif x > 0 && x <= 1
        P4 = x - x^3/factorial(3);
        y = x^(-1/4) * (sin(x) - P4);
    else
        y = 0;
    end
end