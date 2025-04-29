function plotIRFs(IRFMat, m, irf_horizon, shock, shock_horizon, series)
            % Plot IRFs from IRFMat (nx*T×nexog*T)
            if nargin<3
                irf_horizon=80; 
            end
            if nargin<4 
                shock=1;
            end
            if nargin<5
                shock_horizon=0; 
            end
            if nargin<6
                series=m.vars.endo; 
            end
            % determine shock index
            if isnumeric(shock)
                sidx=shock;
            else
                exn=m.vars.exog; 
                sidx=find(strcmp(exn,shock)); 
            end
            vec = IRFMat(:,(sidx-1)*m.T + shock_horizon+1);
            blocks = wide(vec,m.nx);
            figure;
            for i=1:m.nx
                subplot(ceil(m.nx/3),3,i);
                plot(1:irf_horizon, blocks{i}(1:irf_horizon),'LineWidth', 2);
                xlabel('时间');
                title(series{i});
                sgtitle('脉冲响应函数');
            end
        end