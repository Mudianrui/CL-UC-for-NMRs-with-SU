%PPC_rho_inf_derivative
function [y] = PPC_rho_inf_derivative(t,T,epsilon,n,flag)
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
            y = p*(T/t-1)^(p-1)*(-T*t^(-2));
        else
            y = 0;
        end
    otherwise
        error(['Unhandled flag from PPC_rho_inf_derivative = ',num2str(flag)]);
end
