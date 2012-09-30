clear all;
clf;

addpath('test_functions');
adm = @IsAdmissible;% @ ia function handler
nstart = 50;
nend = 65;
maxiters = -1;
maxsize = 8 ^ 2;

x = [];
y = [];

for i = nstart : nend
    A = delsq(numgrid('S', i));
    fprintf('%d %d iterating...\t', i, size(A,1));
    SM = supermatrix(A);
    SM = SM.fulliterate(adm, maxiters, maxsize, 10^-3);
    fprintf('DONE calculating...\t');
    [nnzrk, nodes, total, xnew, ynew] = calc_rkeigenvalues_rec(SM);
    x = [x xnew]; %#ok<AGROW>
    y = [y ynew]; %#ok<AGROW>
    fprintf('DONE %d rk nodes on a total of %d\n', nodes, total);
end

fprintf('making data mmonotonic\t');
[xsorted, ysorted, xx, yy] = make_monotonic(x, y);
fprintf('OK!\n');

fprintf('Interpolating\t');
xxi = 0:10^-1:1;
yyi = interp1([0 xx], [1 yy], xxi, 'spline');
fprintf('OK!\n');
hold on;
plot(xx, yy, 'o');
p = plot(xxi, yyi);
hold off;
set(p, 'Color','red', 'LineWidth',2)
str = sprintf('%d nnz (%d rk nodes) on a total of %d\n', nnzrk, nodes, total);
fprintf(str);
title(str);
