function [e]=find_max_mixture(bo_o,alpha_o,bo_i,alpha_i,start)

lambda_o = 1/sqrt(2*alpha_o);
lambda_i = 1/sqrt(2*alpha_i);

k=1;
e(k) = start;
maxnumsteps = 10000;
while(k<maxnumsteps)
    if(e(k)>lambda_o && e(k)<1-lambda_i)
        if(start==0)
            e(k+1) = lambda_o - 1e-3;
        else
            e(k+1) = 1-lambda_i + 1e-3;
        end
    else
        xi_d_o  = -2*alpha_o * bo_o * e(k) * exp(-alpha_o*e(k)^2);                                  % Eq. (64)
        xi_dd_o =  2*alpha_o * bo_o * ( 2*alpha_o*e(k)^2-1 ) * exp(-2*alpha_o*e(k)^2);              % Eq. (65)
        xi_d_i  =  2*alpha_i * bo_i * (1-e(k)) * exp(-alpha_i*(1-e(k))^2);                          % Eq. (66)
        xi_dd_i =  2*alpha_i * bo_i * ( 2*alpha_i*(1-e(k))^2-1 ) * exp(-2*alpha_i*(1-e(k))^2);      % Eq. (67)
        stepsize = -(xi_d_o+xi_d_i)/(min([0 xi_dd_o]) + min([0 xi_dd_i]));                          % Eq. (68)
        if(abs(stepsize)>50)
            disp('large step')
        end
        if( abs(stepsize)<1e-10 )
            break;
        end
        e(k+1) = e(k) + stepsize;
    end
    k=k+1;
    
end

