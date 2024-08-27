%DirectPPC_amplifier_derivative
function [mua,upsilona] = DirectPPC_amplifier_derivative(z,rho,drho,delta,alpha,n,flag)
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
    case 0
        p=ceil(n+1);
        eta = (rho^2-z^2)/delta^2;
        if z^2<=rho^2-delta^2
            gamma = 1;
        else
            gamma = 1/(1-(eta-1)^p);
        end
    case 4
        eta = ((z/rho)^2-delta^2)/(1-delta^2);
        Deta_z = 2*z/(rho^2*(1-delta^2));
        Deta_rho = -2*z^2/(rho^3*(1-delta^2));
    otherwise
        error(['Unhandled flag from DirectPPC_amplifier = ',num2str(flag),' ',num2str(z),' ',num2str(rho)]);
end

if flag>=1 && flag<=4
    if eta<=0
        gamma = 1;
        Dgamma_eta = 0;
    else
        gamma = 1/(1-eta^(n+1));
        Dgamma_eta = (n+1)*eta^n*gamma^2*delta;
    end
end
Dgamma_z = Dgamma_eta*Deta_z;
Dgamma_rho = Dgamma_eta*Deta_rho;
mu = Dgamma_z+gamma;
upsilon = Dgamma_rho*z*drho;
if eta<=0
    mu = 1;
    upsilon=0;
end
mua = alpha*(mu-1)+1;
upsilona = alpha*upsilon;
