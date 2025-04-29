 function [Xss, Ess] = longsteadystate(m)
            % Long-format steady state: repeat init/terminal and exog across T
            Xss = repmat(m.ss.terminal(:), m.T, 1);
            Ess = repmat(m.ss.exog(:), m.T, 1);
        end