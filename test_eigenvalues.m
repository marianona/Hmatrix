clear all;
hold all;

nstart = 10;
nend = 50;

for n = nstart : nend
    A = delsq(numgrid('S', n));
    ss = svd(full(A));
    fprintf('%d %d %d\n', n, size(A,1), length(ss));
    % normalize
    ss = ss / ss(1);
    lenss = length(ss);
    plot((1:lenss)/lenss, ss);
end
hold off;