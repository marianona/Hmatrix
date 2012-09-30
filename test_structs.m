clear all;
addpath('structs');

A = delsq(numgrid('S', 10));
Afull = full(A);

fullstruct = fullmatrix_struct(A);
rkstruct = rkmatrix_struct(A);
sustruct = supermatrix_struct(A);

full_matrix = fullmatrix(A);
rk_matrix = rkmatrix(A);
su_matrix = supermatrix(A);
whos
