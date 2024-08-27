%PPC_rho_inf
function [y] = PPC_rho_inf(t,T,epsilon,n,flag)
if t<0
    flag = -1;
elseif T<0
    flag = -2;
elseif epsilon<0
    flag = -3;
elseif n~=fix(n) || n<1
    flag = -4;
end

switch flag
    case 0
        p = n+1;
        if t<=T
            y = (T/t-1)^p+epsilon;
        else
            y = epsilon;
        end
    otherwise
        error(['Unhandled flag from PPC_rho_inf = ',num2str(flag)]);
end
