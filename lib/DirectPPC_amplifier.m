%DirectPPC_amplifier
function [gamma] = DirectPPC_amplifier(z,rho,delta,alpha,n,flag)
if abs(z)>=abs(rho)
    flag = -1;
elseif rho<0
    flag = -2;
elseif delta<0
     flag = -3;
elseif delta>0
    if flag==4
        if delta>1
            flag = -3;
        end
    else
        if delta>rho
            flag = -3;
        end
    end
elseif alpha<0
    flag = -4;
elseif n~=fix(n) || n<1
    flag = -5;
end

switch flag
    case 4
        eta = ((z/rho)^2-delta^2)/(1-delta^2);
    otherwise
        error(['Unhandled flag from DirectPPC_amplifier = ',num2str(flag),' ',num2str(z),' ',num2str(rho)]);
end

if flag>=1 && flag<=4
    if eta<=0
        gamma = 1;
    else
        gamma = 1/(1-eta^(n+1));
    end
end
gamma = alpha*(gamma-1)+1;
