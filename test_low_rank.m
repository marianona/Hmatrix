adm = @IsAdmissible;

for i = 30:50
    A = delsq(numgrid('S', i));
    s = supermatrix(A);
    s = s.fulliterate(adm, -1, 256^2, 10^-3);
    fprintf('i: %d\tsize: %d\tnorm: %e\n', i, size(A,1), norm(A - s.getTable));
end