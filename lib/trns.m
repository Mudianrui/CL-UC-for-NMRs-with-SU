%trns
function [s] = trns(x,x1,x2,flag)
if x1==x2
    flag = -1;
    error(['Unhandled flag from trns = ',num2str(flag)]);
end

xx = (x-x1)/(x2-x1);
if xx<=0
    s=0;
elseif xx>=1
    s=1;
else
    switch flag
        case 0
            s = 1/(exp((1-2*xx)/(xx*(1-xx)))+1);
        otherwise
            error(['Unhandled flag from trns = ',num2str(flag)]);
    end
end
