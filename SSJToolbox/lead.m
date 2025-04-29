function v_f = lead(v)
            v_f        = circshift(v,-1);
            v_f(end) =v(end);      % 末行置稳态
        end