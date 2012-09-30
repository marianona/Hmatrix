function result = addparts2_rkmatrix( r, rk_no, rk_mo, rk_r )
    [rows cols] = r.getsize();
    [rk_rrows rk_rcols] = size(rk_r);
    for i = 1:rk_rrows
        for j = 1:rk_rcols
            crnt = rk_r(i,j);
            crntkt = size(crnt.a,2);
            before = rk_no(i,j);
            after = rows - crnt.rows - before;
            crnt.a = [  zeros(before, crntkt)
                        crnt.a
                        zeros(after, crntkt)];
            
            before = rk_mo(i,j);
            after = cols - crnt.cols - before;
            crnt.b = [  zeros(before, crntkt)
                        crnt.b
                        zeros(after, crntkt)];
            crnt.rows = rows;
            crnt.cols = cols;
            rk_r(i,j) = crnt;
        end
    end
    for i = 1:numel(rk_r)
        r = add_rkmatrix(r, r, rk_r(i));
    end
    result = r;
end
