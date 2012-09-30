clear all;
clf;
addpath('test_functions');
adm = @IsAdmissible;% @ ia function handler
useFull = 1;
nstart = 20;
nend = 50;
maxiters = -1;
maxsize = 8 ^ 2;
relativeError = 0;
per(nend) = 0;
eigenper(nend) = 0;
memN(nend) = 0;
memH(nend) = 0;

for i = nstart : nend
    if useFull == 1
        A = full(delsq(numgrid('S', i)));
    else
        A = delsq(numgrid('S', i));
    end
    wh = whos('A');
    memN(i) = wh.bytes;
    fprintf('%d %d iterating...\t', i, size(A,1));
    SM = supermatrix(A);
    SM = SM.fulliterate(adm, maxiters, maxsize, relativeError);
    wh = whos('SM');
    memH(i) = wh.bytes;
    fprintf('DONE calculating...\t');
    [total, rk, rkper] = stat_rk(SM);
    per(i) = rk / total;
    eigenper(i) = rkper / rk;
    fprintf('DONE %d rk nodes on a total of %d with %f percantage of nonzero eigen values\n', rk, total, rkper / rk);
    clear A SM;
end

subplot(1,2,1);
plot(nstart:nend, per(nstart:nend),...
     nstart:nend, eigenper(nstart:nend));
legend('rk nodes percentage','non-zero eigen values percentage');
xlabel('delsq numgrid parameter');

subplot(1,2,2);
plot(nstart:nend, memN(nstart:nend),...
     nstart:nend, memH(nstart:nend));
legend('Matlab','Hmatrix');
xlabel('delsq numgrid parameter');
ylabel('memory');
title('Memory Usage Comparison');
