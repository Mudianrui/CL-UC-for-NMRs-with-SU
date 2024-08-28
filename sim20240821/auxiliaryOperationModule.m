function [sys,x0,str,ts] = auxiliaryOperationModule(t,x,u,flag)
switch flag
case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
case 3
    sys=mdlOutputs(t,x,u);
case {1,2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 5;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [];

function sys=mdlOutputs(~,~,u)
xd=u(1);
yd=u(2);
thd=u(3);
x=u(4);
y=u(5);
th=u(6);

eth=thNormalization(th-thd);

global Delta_x Delta_y
global lc
%Control point coordinates
Psi = [cos(th),-sin(th);sin(th),cos(th)];
ps=[x;y];
p0=ps-Psi*[Delta_x;Delta_y];
xc=p0(1)+lc*cos(th);
yc=p0(2)+lc*sin(th);
%Actual error of control points
xeT=xc-xd;
yeT=yc-yd;

sys(1)=eth;
sys(2)=xc;
sys(3)=yc;
sys(4)=xeT;
sys(5)=yeT;


