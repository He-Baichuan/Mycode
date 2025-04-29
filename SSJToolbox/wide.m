 function blocks = wide(x, n)
            % Split long vector x (n*T×1) into cell array of n vectors length T
            T = numel(x)/n;
            x = reshape(x,n,T);
            x = x';
            blocks = mat2cell(x,T, ones(1,n));
        end