Hmatrix
=======

A simple library for Matlab, for experementing with Hierarchical Matrices theory.

The implementation is based on the [paper](http://www.mis.mpg.de/de/publications/andere-reihen/ln/lecturenote-2103.html) written by Steffen Börm, Lars Grasedyck, and Wolfgang Hackbusch.

##Simple demonstration code:
```matlab
% This should be the path to the library.
% Let as is if you placed the Hmatrix folder next to your *.m files.
addpath('Hmatrix');

A = delsq(numgrid('S', 10)); % assing YOUR table here (should show A = mytable; instead).

adm = @IsAdmissible; % this points to the default admissibility condition.
maxiterations = -1;
minBlockSize = 256^2;
relativeError = 0;

% initializes a Hmatrix tree with depth 1.
S = supermatrix(full(A));
% Does the actual tree structuring.
S = S.fulliterate(adm, maxiterations, minBlockSize, relativeError);

result = S.invert(); % result holds the inverse of table A.
```

More simple examples, can be found in main.m file.