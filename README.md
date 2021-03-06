Hmatrix
=======

A simple library for Matlab, for experementing with Hierarchical Matrices theory.

The implementation is based on the [paper](http://www.mis.mpg.de/de/publications/andere-reihen/ln/lecturenote-2103.html) written by Steffen Börm, Lars Grasedyck, and Wolfgang Hackbusch.

##Simple demonstration code:
```matlab
addpath('Hmatrix');% This should be the path to the library.
% Let as is if you placed the Hmatrix folder next to your *.m files.

A = delsq(numgrid('S', 10)); % assing YOUR table here (should show A = mytable; instead).

adm = @IsAdmissible; % this points to the default admissibility condition.
maxiterations = -1;
minBlockSize = 256^2;
relativeError = 0;

S = supermatrix(full(A)); % initializes a Hmatrix tree with depth 1.
S = S.fulliterate(adm, maxiterations, minBlockSize, relativeError); % Does the actual tree structuring.

result = S.invert(); % result holds the inverse of table A in supermatrix form.
result.getTable() % returns a simple matlab matrix from the supermatrix representation.

```

More simple examples, can be found in main.m file.