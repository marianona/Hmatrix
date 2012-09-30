function r = addprod2_rkmatrix( r, a, b )
    if (a.cols ~= b.rows || r.rows ~= a.rows || r.cols ~= b.cols)
        error('Matrix dimensions must agree.\n');
    end
    bn = a.block_rows;
    bam = a.block_cols;
    bm = b.block_cols;
    if a.hasS && b.hasS
        rk_no = zeros(bn, bm);
        rk_mo = zeros(bn, bm);
        rk_r(bn, bm) = rkmatrix;
        no = 0;
        for i = 1:bn
            mo = 0;
            for j = 1:bm
                rk_no(i, j) = no;
                rk_mo(i, j) = mo;
                tmp = rkmatrix(zeros(a.s(i,1).rows,b.s(1,j).cols));
                tmp.k = r.k;
                rk_r(i, j) = tmp;
                for k = 1:bam
                    rk_r(i, j) = addprod2_rkmatrix(rk_r(i, j), a.s(i, k), b.s(k, j));
                end
                mo = mo + b.s(1,j).cols;
            end
            no = no + a.s(i,1).rows;
        end
        r = addparts2_rkmatrix(r, rk_no, rk_mo, rk_r);
    else
        if a.hasR && b.hasR
            tmp = a.r * b.r;
        else
            tmp = rkmatrix(a.getTable * b.getTable);
        end
        r = add_rkmatrix(r, r, tmp);
    end
end
