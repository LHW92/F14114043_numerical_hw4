clc; clear; close all;

% 定義函數
f = @(x) (x.^2) .* log(x);

% 設定積分區間
a = 1;
b = 1.5;
h = (b - a) / 2; % 變數變換的係數

% Gauss-Legendre nodes and weights for n=3
t_3 = [-sqrt(3/5), 0, sqrt(3/5)];
c_3 = [5/9, 8/9, 5/9];

% Gauss-Legendre nodes and weights for n=4
t_4 = [-0.861136, -0.339981, 0.339981, 0.861136];
c_4 = [0.347855, 0.652145, 0.652145, 0.347855];

% 轉換節點
x_3 = h * t_3 + (a + b) / 2;
x_4 = h * t_4 + (a + b) / 2;

% 計算 f(x) 值
fx_3 = f(x_3);
fx_4 = f(x_4);

% 計算高斯積分
I_3 = h * sum(c_3 .* fx_3);
I_4 = h * sum(c_4 .* fx_4);

syms x;
I_exact = int(x^2 * log(x), 1, 1.5);
vpa(I_exact, 6);  % 取6位小數

% Compare these results with the exact value
e3_abs = abs(I_3 - I_exact);
e3_rel = abs(e3_abs/I_exact);
e4_abs = abs(I_4 - I_exact);
e4_rel = abs(e4_abs/I_exact);

% 顯示結果
fprintf('Exact Integral Value: %.6f\n', I_exact);
fprintf('---------------------------------------------------\n');
fprintf('Gaussian Quadrature with n = 3: %.6f\n', I_3);
fprintf("Absolute Error = %.6f, Relative Error = %.4f%%\n", e3_abs, e3_rel*100);
fprintf('---------------------------------------------------\n');
fprintf('Gaussian Quadrature with n = 4: %.6f\n', I_4);
fprintf("Absolute Error = %.6f, Relative Error = %.4f%%\n", e4_abs, e4_rel*100);