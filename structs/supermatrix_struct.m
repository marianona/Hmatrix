function [ output_args ] = supermatrix_struct( table )
    if numel(table) == 1
        [r,c] = size(table);
        br = 0;
        bc = 0;
    else
        br = size(table, 1);
        bc = size(table, 2);
        r = 0;
        c = 0;
        for i = 1:br
            r = r + size(table(i,1), 1);
        end
        for i = 1:bc
            c = c + size(table(1,i), 2);
        end
    end
    
    output_args = struct( ...
        'rows', r, ...
        'cols', c, ...
        'block_rows', br, ...
        'block_cols', bc ...
    );
    if numel(table) == 1
        output_args.f = table;
    else
        output_args.s = table;
    end
end
