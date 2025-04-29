function v_l = lag(v)
            v_l      = circshift(v,1);
            v_l(1) = v(1);        % 第一行置稳态
        end