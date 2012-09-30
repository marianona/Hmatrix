function c = mul_supermatrix( a, b )
    if (a.cols ~= b.rows )
        error('Matrix dimensions must agree.\n');
    end
    bn = a.block_rows;
    bm = b.block_cols;
    bam = a.block_cols;
    if a.hasS && b.hasS
        c = supermatrix();
        c.rows = a.rows;
        c.cols = b.cols;
        c.block_rows = bn;
        c.block_cols = bm;
        c.s = supermatrix;
        c.s(bn, bm) = supermatrix;
        for i = 1:bn
            for j = 1:bm
                for k = 1:bam
                    if c.s(i,j).rows == 0
                        c.s(i,j) = mul_supermatrix(a.s(i,k), b.s(k,j));
                    else
                        c.s(i,j) = c.s(i,j) + mul_supermatrix(a.s(i,k), b.s(k,j));
                        %c.s(i,j) = muladd_supermatrix(c.s(i,j), a.s(i,k), b.s(k,j));
                    end
                end
            end
        end
    elseif a.hasR && b.hasR
        tmp = a.r * b.r;
        c = tmp.getsupermatrix();
    elseif a.hasF && b.hasF
        tmp = a.f * b.f;
        c = tmp.getsupermatrix();
    else
        c = supermatrix(a.getTable * b.getTable);
    end
end
