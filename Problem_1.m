clc; clear; close all;

% 定義被積函數
f = @(x) exp(x) .* sin(4*x);

% 積分區間
a = 1;
b = 2;

% 分割區間數
n = 5;
h = (b - a) / (2*n);
x = linspace(a, b, 2*n+1);

% 1. Composite Trapezoidal Rule
trap_result = (h/2) * (f(x(1)) + 2*sum(f(x(2:end-1))) + f(x(end)));

% 2. Composite Simpson's Rule
if mod(n,2) == 1
    n = n + 1; % 辛普森法則需要偶數個子區間
end
simpson_result = (h/3) * (f(x(1)) + 4*sum(f(x(2:2:end-1))) + 2*sum(f(x(3:2:end-2))) + f(x(end)));

% 3. Composite Midpoint Rule
midpoints = a+h:2*h:b-h;
midpoint_result = 2 * h * sum(f(midpoints));

% 解析解（作為誤差比較）
syms x;
I_exact = int(exp(x) * sin(4*x), x, a, b);
I_exact_value = double(I_exact);

% 顯示結果
fprintf('Exact Integral Value: %.6f\n', I_exact_value);
fprintf('a. Composite Trapezoidal Rule Result: %.6f\n', trap_result);
fprintf('b. Composite Simpson''s Rule Result: %.6f\n', simpson_result);
fprintf('c. Composite Midpoint Rule Result: %.6f\n', midpoint_result);