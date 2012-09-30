%#ok<*NOPTS>
adm = @IsAdmissible;
A = delsq(numgrid('S', 6));
s = supermatrix(A);
s = s.iterate(adm);
s = s.iterate(adm);
maxerr = 10^-6;

fprintf('rsvd\n');
B = [1 2 3 4
    5 6 7 8
    9 10 11 12];
r1 = rkmatrix(B);
[u1 s1 v1]=r1.rsvd();
test1 = u1*s1*v1';
if norm(B - test1) <= maxerr
    fprintf('OK!\t');
else
    test1
end

r2 = rkmatrix(B');
[u1 s1 v1]=r2.rsvd();
test2 = u1*s1*v1';
if norm(B' - test2) <= maxerr
    fprintf('OK!\t');
else
    test2
end

C = [1 2 0
    3 4 0
    0 0 0];
r3 = rkmatrix(C);
[u1 s1 v1]=r3.rsvd();
test3 = u1*s1*v1';
if norm(C - test3) <= maxerr
    fprintf('OK!\t');
else
    test3
end

D = [0 0 0
    0 1 2
    0 3 4];
r4 = rkmatrix(D);
[u1 s1 v1]=r4.rsvd();
test4 = u1*s1*v1';
if norm(D - test4) <= maxerr
    fprintf('OK!\t');
else
    test4
end

rzero = rkmatrix(zeros(3));
[u1 s1 v1]=rzero.rsvd();
test5 = u1*s1*v1';
if norm(test5) <= maxerr
    fprintf('OK!\t');
else
    test5
end
fprintf('\n');

fprintf('add_rkmatrix r3 r4\n');
result = r3.getTable() + r4.getTable();
r34=add_rkmatrix(r3,r3,r4);
test1 = r34.getTable;
if norm(result - test1) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test1
end
r34o=r3+r4;
test2 = r34o.getTable;
if norm(result - test2) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test2
end

fprintf('\n');

fprintf('add_rkmatrix r3 rzero\n');
result = r3.getTable() + rzero.getTable();
r35=add_rkmatrix(r3,r3,rzero);
test1 = r35.getTable;
if norm(result - test1) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test1
end
r35o=r3+rzero;
test2 = r35o.getTable;
if norm(result - test2) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test2
end
fprintf('\n');

fprintf('add_rkmatrix rzero rzero\n');
rzero = add_rkmatrix(rzero,rzero,rzero);
test1 = rzero.getTable;
if norm(test1) <= maxerr
    fprintf('OK!\t');
else
    test1
end
fprintf('\n');

fprintf('add_rkmatrix a b\n');
a = rkmatrix(s.s(1).getTable());
b = rkmatrix(s.s(2).getTable());
result = a.getTable + b.getTable;
%full(result)
c0 = rkmatrix();
c0.k = 8;
c = add_rkmatrix(c0,a,b);
test1 = c.getTable();
if norm(result - test1) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test1
end
c=a+b;
test2 = c.getTable();
if norm(result - test2) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test2
end
fprintf('\n');

a = supermatrix(s.s(1).getTable());
a = a.iterate(adm);
a = a.iterate(adm);
b = supermatrix(s.s(2).getTable());
b = b.iterate(adm);
b = b.iterate(adm);

fprintf('supermatrix plus a b\n');
result = a.getTable + b.getTable;
%full(result)
d = a + b;
test1 = full(d.getTable);
if norm(result - test1) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test1
end
fprintf('\n');

fprintf('addprod2_rkmatrix a b\n');
result = a.getTable * b.getTable;
%full(result)
d0 = rkmatrix(zeros(a.rows, b.cols));
d0.k = 8;
d = addprod2_rkmatrix(d0, a, b);
test1 = d.getTable;
if norm(result - test1) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test1
end
fprintf('\n');

fprintf('muladd_supermatrix a b\n');
result = a.getTable * b.getTable;
%full(result)
d0 = supermatrix(zeros(a.rows, b.cols));
d = muladd_supermatrix(d0, a, b);
test1 = full(d.getTable);
if norm(result - test1) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test1
end
d = a * b;
test2 = full(d.getTable);
if norm(result - test2) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    test2
end
fprintf('\n');

fprintf('invert\n');
d=a;
result = inv(d.getTable);
%full(result)
result2 = d.invert;
test1 = full(result2.getTable);
if norm(result - test1) <= maxerr
    fprintf('OK!\t');
else
    full(result)
    full(test1)
    fprintf('invert rounded diff\n');
    round(full(result - result2.getTable)*10^10)/10^10
end
fprintf('\n');
