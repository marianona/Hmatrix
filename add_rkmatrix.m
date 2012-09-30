function c = add_rkmatrix(c, a, b)
    if (a.rows ~= b.rows || a.cols ~= b.cols)
        error('Matrix dimensions must agree.\n');
    end
    kt = a.kt + b.kt;
    a_plus_b = rkmatrix();
    a_plus_b.k = kt;
    a_plus_b.kt = kt;
    [a_plus_b.rows a_plus_b.cols] = a.getsize();
    a_plus_b.a = [a.a b.a];
    a_plus_b.b = [a.b b.b];
    
    [u s v] = a_plus_b.rsvd();
    s = diag(s);% does not contain zeros, so length == nnz

    if isa(c.tc, 'truncation_control')% && c.tc.adaptive == 1
        maxiters = min(kt, length(s));
        i = 0;
        while i < maxiters && s(i+1) > c.tc.abs_eps && s(i+1) > c.tc.rel_eps * s(1)
            i = i + 1;
        end
        if i > c.k
            c.k = i;
        end
    else
        i = min(kt, c.k);
    end
    i = min(i, length(s));
    c.kt = i;
    [c.rows c.cols] = a.getsize();
    c.a = u(:, 1:i);
    %c.b = v(:, 1:i) .* (ones(a.cols, 1) * s(1:i)');
    c.b = v(:, 1:i) * diag(s(1:i));
end
