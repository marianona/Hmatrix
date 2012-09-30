function [ total, numOfRK, rkeigen_perc ] = stat_rk( SM )
    total = 0;
    numOfRK = 0;
    rkeigen_perc = 0;
    if SM.hasR
        ss = svd(full(SM.getTable));
        n = 0;
        while n+1 <= length(ss) && ss(n+1) ~= 0
            n = n + 1;
        end
        numOfRK = 1;
        total = 1;
        rkeigen_perc = n / length(ss);
    elseif SM.hasS
        for i = 1 : numel(SM.s)
            [tot, rk, rkper] = stat_rk(SM.s(i));
            total = total + tot;
            numOfRK = numOfRK + rk;
            rkeigen_perc = rkeigen_perc + rkper;
        end
    else
        total = 1;
    end
end

