% this script depends on UFget matlab library
clc;
clear all;
clf;
adm = @IsAdmissible;
useFull = 1;
testloops = 1;
minBlockSize = 256^2;
relativeError = 0;
showerror = 0;
stopAfterNExceptions = 3; % stop after # in row exceptions.

%% Matrices selected from UFgui:
% Example usage:
% for k = 1:length(ids)
%    Problem = UFget (ids (k))
% end
ids = [
13 % HB/bcspwr01
14 % HB/bcspwr02
15 % HB/bcspwr03
16 % HB/bcspwr04
17 % HB/bcspwr05
2 % HB/494_bus
3 % HB/662_bus
4 % HB/685_bus
1627 % Bai/qh768
328 % Bai/qh882
1 % HB/1138_bus
140 % HB/eris1176
18 % HB/bcspwr06
1626 % Bai/qh1484
19 % HB/bcspwr07
20 % HB/bcspwr08
21 % HB/bcspwr09
2230 % TSOPF/TSOPF_FS_b9_c1
155 % HB/gemat1
156 % HB/gemat11
157 % HB/gemat12
22 % HB/bcspwr10
2232 % TSOPF/TSOPF_RS_b162_c1
2245 % TSOPF/TSOPF_RS_b9_c6
2539 % IPSO/TSC_OPF_1047
2540 % IPSO/TSC_OPF_300
2221 % TSOPF/TSOPF_FS_b162_c1
2242 % TSOPF/TSOPF_RS_b39_c7
1871 % YCheng/psse1
2214 % QY/case9
2231 % TSOPF/TSOPF_FS_b9_c6
2237 % TSOPF/TSOPF_RS_b300_c1
2233 % TSOPF/TSOPF_RS_b162_c3
2537 % IPSO/OPF_3754
1212 % LiuWenzhuo/powersim
2243 % TSOPF/TSOPF_RS_b678_c1
2234 % TSOPF/TSOPF_RS_b162_c4
1874 % HVDC/hvdc1
2235 % TSOPF/TSOPF_RS_b2052_c1
1870 % YCheng/psse0
2229 % TSOPF/TSOPF_FS_b39_c7
2238 % TSOPF/TSOPF_RS_b300_c2
1872 % YCheng/psse2
2220 % TSOPF/TSOPF_FS_b300
2224 % TSOPF/TSOPF_FS_b300_c1
2538 % IPSO/OPF_6000
2222 % TSOPF/TSOPF_FS_b162_c3
2244 % TSOPF/TSOPF_RS_b678_c2
2240 % TSOPF/TSOPF_RS_b39_c19
2219 % TSOPF/TSOPF_RS_b2383
2236 % TSOPF/TSOPF_RS_b2383_c1
2136 % QY/case39
2223 % TSOPF/TSOPF_FS_b162_c4
2239 % TSOPF/TSOPF_RS_b300_c3
2536 % IPSO/OPF_10000
2225 % TSOPF/TSOPF_FS_b300_c2
2241 % TSOPF/TSOPF_RS_b39_c30
2227 % TSOPF/TSOPF_FS_b39_c19
2226 % TSOPF/TSOPF_FS_b300_c3
2228 % TSOPF/TSOPF_FS_b39_c30
1875 % HVDC/hvdc2
2534 % IPSO/HTC_336_4438
2535 % IPSO/HTC_336_9129
] ;

nstart = 2;
nend = length(ids);
        
t1(nend) = 0;
t2(nend) = 0;
tbuild(nend) = 0;
err(nend) = 0;
catches(nend) = 0;

i = nstart;
nExceptios = 0;
while true
    if i > nend
        i = i - 1;
        break;
    end
    try
        Problem = UFget(ids (i));
        if useFull == 1
            A = full(Problem.A);
        else
            A = Problem.A;
        end
        
        s = size(A);
        if s(1) * s(2) < minBlockSize
            i = i + 1;
            nstart = i;
            continue;
        end
        
        fprintf('i: %d/%d size: %d\t', i, length(ids), size(A,1));
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
        t1(i) = tt1;
        fprintf('normal: %f\t', tt1);
        
        ttbuild = 0;
        tic;
        s = supermatrix(A);
        clear A;
        s.fulliterate(adm, -1, minBlockSize, relativeError);
        ttbuild = toc;
        tbuild(i) = ttbuild;
        
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
        t2(i) = tt2;
        fprintf('Hmatrix: %f\t', tt2);
        fprintf('build time: %f\t', ttbuild);
        clear tt1 tt2 ttbuild;
        
        if showerror == 1
            if issparse(res1)
                err(i) = norm(full(res2.getTable - res1)) / norm(full(res1));
            else
                err(i) = norm(res2.getTable - res1) / norm(res1);
            end
            fprintf('error: %e', err(i));
            clear res1 res2;
         else
             err(i) = 0;
         end
        
        fprintf('\n');
        nExceptions = 0; % reset the exception counter.
    catch ex
        fprintf('\n!!!Exception!!!\n');
        nExceptions = nExceptions + 1;
        if nExceptions >= stopAfterNExceptions
            i = i - nExceptions;
            break;
        end
    end
    i = i + 1;
end

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
