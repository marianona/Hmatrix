function c = muladd_supermatrix( c, a, b )
    if (a.cols ~= b.rows || c.rows ~= a.rows || c.cols ~= b.cols)
        error('Matrix dimensions must agree.\n');
    end
    bn = c.block_rows;
    bm = c.block_cols;
    bam = a.block_cols;
    if a.hasS && b.hasS
        if c.hasR
            c.r = addprod2_rkmatrix(c,r, a, b);
        elseif c.hasF
            tmp = mul_supermatrix(a, b);
            c.f.e = c.f.e + tmp.getTable;
        elseif c.hasS
            for i = 1:bn
                for j = 1:bm
                    for k = 1:bam
                        c.s(i,j) = muladd_supermatrix(c.s(i,j), a.s(i,k), b.s(k,j));
                    end
                end
            end
        end
    elseif a.hasR && b.hasR && c.hasR
        %c.r = c.r + a.r * b.r;
        %c.r.tc = truncation_control();
        c.r = add_rkmatrix(c.r, c.r, a.r * b.r);
    elseif (a.hasR && b.hasR) || (a.hasF && b.hasF)
        tmp = a * b;
        c = supermatrix(c.getTable + tmp.getTable);
    else
        c = supermatrix(c.getTable + a.getTable * b.getTable);
    end
end
