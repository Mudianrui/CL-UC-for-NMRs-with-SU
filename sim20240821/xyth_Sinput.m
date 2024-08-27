function [sys,x0,str,ts] = xyth_Sinput(t,x,u,flag)
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
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [];

addpath(genpath('..\lib'))
%controllerTypeK：  52：Traditional adaptive controller（Comparison）
%                   54：Unbiased Controller
%                   56：Unbiased controller with prescribed performance
global controllerTypeK
controllerTypeK = 52;
global Delta_x Delta_y
Delta_x = 0.5;
Delta_y = 0.1;
global lc
lc = 0.1;

function sys=mdlOutputs(~,~,u)
t = u(1);
wd=0.5;
Ax=1;
Ay=1;
N=4;
sys(1)=Ay*(cos(wd*t)+0.5*cos(N*wd*t));%x
sys(2)=Ax*(sin(wd*t)+0.5*sin(N*wd*t));%y
dx = -Ax*(wd*sin(wd*t)+0.5*N*wd*sin(N*wd*t));
dy = Ax*(wd*cos(wd*t)+0.5*N*wd*cos(N*wd*t));
sys(3)=atan2(dy,dx);%th
sys(4)=dx;%dx
sys(5)=dy;%dy
