function S = id(name,m,num)
            if num == 0
                S = m.idx_now(name);
            
            end
            if num == -1
                S= m.idx_now([name '_l']);
               
            end
            if num ==1
                S = m.idx_now([name '_p']);
               
            end
        end