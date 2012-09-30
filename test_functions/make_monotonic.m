function [ x, y, xx, yy ] = make_monotonic( x, y )
    % selection sort
    for i = 1:(numel(x)-1)
        fprintf('sorting: %f\n', 100 * i / (numel(x)-1));
        min = x(i);
        minpos = i;
        for j = (i+1):numel(x)
            if x(j) < min
               min = x(j); 
               minpos = j;
            end
        end
        if min < x(i)
            tmp = x(i);
            x(i) = x(minpos);
            x(minpos) = tmp;
            
            tmp = y(i);
            y(i) = y(minpos);
            y(minpos) = tmp;
        end
    end
    
    % average y values on the same x
    xx = [];
    yy = [];
    i = 1;
    while i <= numel(x)
        fprintf('averaging: %f\n', 100 * i / numel(x));
        prev = x(i);
        j = 0;
        while i + j + 1 <= numel(x) &&  prev == x(i + j + 1)
            j = j + 1;
        end
        avg = sum(y(i:(i+j))) / (j+1);
        xx = [xx x(i)]; %#ok<AGROW>
        yy = [yy avg];  %#ok<AGROW>
        i = i + j + 1;
    end
    
%     % fix x duplicates by adding eps
%     for i = 1:numel(x)
%         fprintf('de-duplicating: %f\n', 100 * i / numel(x));
%         prev = x(i);
%         j = 0;
%         while i + j + 1 <= numel(x) &&  prev == x(i + j + 1)
%             j = j + 1;
%             fprintf('de-duplicating: %f\n', 100 * (i + j) / numel(x));
%             x(i + j) = x(i + j) + j * eps;
%         end
%         i = i + j; %#ok<FXSET>
%     end
end

