clc;
clear all;
clf;
adm = @IsAdmissible;
useFull = 1;
nstart = 85;
nend = Inf;
testloops = 1;
minBlockSize = 256^2;
relativeError = 0;
showerror = 0;
fprintf('start: %d, end: %d, loops: %d, blocksize: %d, error: %e\n', nstart, nend, testloops, minBlockSize^(1/2), relativeError);

t1(1) = 0;
t2(1) = 0;
tbuild(1) = 0;
err(1) = 0;
i = nstart;
while true
    if i > nend
        break;
    end
    try
        if useFull == 1
            A = full(delsq(numgrid('S', i)));
        else
            A = delsq(numgrid('S', i));
        end
        fprintf('i: %d size: %d\t', i, size(A,1));
        tt1 = 0;
        for ii = 1:testloops
            tic;
            res1 = inv(A);
            tt1 = tt1 + toc;
            if showerror == 0
                clear res1;
            end
        end
        tt1 = tt1 / testloops;
        t1(i) = tt1; %#ok<SAGROW>
        fprintf('normal: %f\t', tt1);
        
        ttbuild = 0;
        tic
        s = supermatrix(A);
        clear A;
        s = s.fulliterate(adm, -1, minBlockSize, relativeError);
        ttbuild = toc;
        tbuild(i) = ttbuild; %#ok<SAGROW>
        
        tt2 = 0;
        for ii = 1:testloops
            tic;
            res2 = s.invert();
            tt2 = tt2 + toc;
            if showerror == 0
                clear res2;
            end
        end
        clear s;
        tt2 = tt2 / testloops;
        t2(i) = tt2; %#ok<SAGROW>
        fprintf('Hmatrix: %f\t', tt2);
        fprintf('build time: %f\t', ttbuild);
        clear tt1 tt2 ttbuild;
        
        if showerror == 1
            if issparse(res1)
                err(i) = norm(full(res2.getTable - res1)) / norm(full(res1)); %#ok<SAGROW>
            else
                err(i) = norm(res2.getTable - res1) / norm(res1); %#ok<SAGROW>
            end
            fprintf('error: %e', err(i));
            clear res1 res2;
         else
             err(i) = 0; %#ok<SAGROW>
         end
        
        fprintf('\n');
    catch ex
        break;
    end
    i = i + 1;
end

i = i - 1;
subplot(2,2,1);
plot(nstart:i, t1(nstart:i),...
     nstart:i, t2(nstart:i),...
     nstart:i, tbuild(nstart:i));
xlabel('delsq numgrid parameter');
ylabel('time in seconds');
title('Time comparison');
legend('normal inversion','Hmatrix inversion','build time');

diff = ((t2(1:i) - t1(1:i)) ./ t1(1:i)) * 100;
subplot(2,2,2);
plot(nstart:i, diff(nstart:i),...
     nstart:i, zeros(1, 1+i-nstart)); % horizontal line
xlabel('delsq numgrid parameter');
ylabel('time difference percentage');
title('Time Difference comparison');

if showerror == 1
    subplot(2,2,3);
    plot(nstart:i, err(nstart:i))
    xlabel('delsq numgrid parameter');
    ylabel('relative error');
    title('Error comparison');
end
