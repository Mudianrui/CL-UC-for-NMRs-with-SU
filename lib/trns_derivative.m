%trns_derivative
function [ds] = trns_derivative(x,x1,x2,flag)
if x1==x2
    flag = -1;
    error(['Unhandled flag from trns = ',num2str(flag)]);
end

xx = (x-x1)/(x2-x1);
if xx<=0 || xx>=1
    ds=0;
else
    switch flag
        case 0
            ds = -1/(exp((1-2*xx)/(xx*(1-xx)))+1)^2 * exp((1-2*xx)/(xx*(1-xx))) * (-2*xx*(1-xx)-(1-2*xx)*(1-2*xx))/(xx^2*(1-xx)^2);
        otherwise
            error(['Unhandled flag from trns = ',num2str(flag)]);
    end
end
