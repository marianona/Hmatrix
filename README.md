Hmatrix
=======

A simple library for Matlab, for experementing with Hierarchical Matrices theory.

The implementation is based on the [paper](http://www.mis.mpg.de/de/publications/andere-reihen/ln/lecturenote-2103.html) written by Steffen Börm, Lars Grasedyck, and Wolfgang Hackbusch.

Simple demonstration code:
```matlab
% this should be the path to the library.
% if your *.m files are in the same folder with the Hmatrix library, there is need for change.
addpath('Hmatrix');

A = delsq(numgrid('S', 10)); % assing YOUR table here

adm = @IsAdmissible; % this points to the default admissibility condition.
maxiterations = -1;
minBlockSize = 256^2;
relativeError = 0;

S = supermatrix(full(A));
S = S.fulliterate(adm, maxiterations, minBlockSize, relativeError);

result = S.invert(); % result holds the inverse of table A.
```