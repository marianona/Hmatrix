function [ output_args ] = rkmatrix_struct( table )
    [r,c] = size(table);
    
    if issparse(table)
        [a s v] = svd(full(table));
    else
        [a s v] = svd(table);
    end
    k = nnz(diag(s));
    if  r > k || c > k
        a = a(:, 1:k);
        b = v(:, 1:k) * s(1:k, 1:k)';
    else
        b = v * s';
    end
    
    output_args = struct( ...
        'k', k, ...
        'kt', k, ...
        'rows', r, ...
        'cols', c, ...
        'a', a, ...
        'b', b ...
    );
end