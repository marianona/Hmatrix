function [ output ] = IsAdmissible( t, s )
    output = max(diam(t), diam(s)) <= dist(t, s);
end

function [ output ] = diam( t )
    output = t(2) - t(1);
    if output < 0
        output = -output;
    end
end

function [ output ] = dist( t, s )
    ts = t(1);
    te = t(2);
    ss = s(1);
    se = s(2);
    if (ts <= ss && ss <= te) || (ts <= se && se <= te) || (ss < ts && te < se)
        output = 0;
    else
        if te < ss
            output = ss - te;
        elseif se < ts
            output = ts - se;
        else
            output = -1;
        end
    end
end