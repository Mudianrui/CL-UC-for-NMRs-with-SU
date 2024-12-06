close all;
load('sim20241206_1_50.mat')
AC_Deltas_hat = Deltas_hat;
AC_th = th;
AC_vdwd = vdwd;
AC_x = x;
AC_y = y;
load('sim20241206_1_52.mat')
PPC_Deltas_hat = Deltas_hat;
PPC_th = th;
PPC_vdwd = vdwd;
PPC_x = x;
PPC_y = y;
load('sim20241206_1_54.mat')
UAC_Deltas_hat = Deltas_hat;
UAC_th = th;
UAC_vdwd = vdwd;
UAC_x = x;
UAC_y = y;
load('sim20241206_1_56.mat')
UPPC_Deltas_hat = Deltas_hat;
UPPC_th = th;
UPPC_vdwd = vdwd;
UPPC_x = x;
UPPC_y = y;

DH = [AC_Deltas_hat,PPC_Deltas_hat,UAC_Deltas_hat,UPPC_Deltas_hat];
X = [AC_x,PPC_x,UAC_x,UPPC_x];
Y = [AC_y,PPC_y,UAC_y,UPPC_y];

FontSize = 12;

%% X-Y
figure('color','white');
grid on
% axes('position',[0.1502 0.2210 0.45 0.7040]);%[0.1502 0.2210 0.7548 0.7040]
hold on
plot(X(1).signals.values(:,1),Y(1).signals.values(:,1),'c-','LineWidth',2.5);
plot(X(1).signals.values(:,6),Y(1).signals.values(:,6),'b--','LineWidth',2.5);
plot(X(2).signals.values(:,6),Y(2).signals.values(:,6),'r--','LineWidth',2);
plot(X(3).signals.values(:,6),Y(3).signals.values(:,6),'k--','LineWidth',1.5);
plot(X(4).signals.values(:,6),Y(4).signals.values(:,6),'g--','LineWidth',1);
scatter(X(1).signals.values(1,6),Y(1).signals.values(1,6),100,'m',"filled")
axis([-1.5 1.8 -1.8 1.5]);
xlabel('$x$(m)','Interpreter','latex','FontSize',FontSize);
ylabel('$y$(m)','Interpreter','latex','FontSize',FontSize);
h=legend('$[x_r,y_r]$','$[x_c,y_c]$(AC in [22])','$[x_c,y_c]$(PPC in [20])','$[x_c,y_c]$(proposed UC)','$[x_c,y_c]$(proposed PPUC)');
set(gca,'FontName','Helvetica','FontSize',FontSize);
set(h,'Interpreter','latex');
% gca
% axes('position',[0.7 0.18 0.25 0.5]);
% hold on
% plot(EOxO(:,1),EOyO(:,1),'r-.','LineWidth',1);
% axis([xe0-re-0.1 xe0-re+0.1 -1.2 1.2]);

%% X
figure('color','white');
grid on
hold on
plot(X(1).time,X(1).signals.values(:,1),'c-','LineWidth',2.5);
plot(X(1).time,X(1).signals.values(:,6),'b-','LineWidth',2);
plot(X(2).time,X(2).signals.values(:,6),'r--','LineWidth',2);
plot(X(3).time,X(3).signals.values(:,6),'k-','LineWidth',2);
plot(X(4).time,X(4).signals.values(:,6),'g--','LineWidth',2);
ylim([-1.5,2]);
xlabel('$t$(s)','Interpreter','latex','FontSize',FontSize);
ylabel('$x$(m)','Interpreter','latex','FontSize',FontSize);
h=legend('$x_r$','$x_c$(AC in [22])','$x_c$(PPC in [20])','$x_c$(proposed UC)','$x_c$(proposed PPUC)');
set(gca,'FontName','Helvetica','FontSize',FontSize);
set(h,'Interpreter','latex');

%% Y
figure('color','white');
grid on
hold on
plot(Y(1).time,Y(1).signals.values(:,1),'c-','LineWidth',2.5);
plot(Y(1).time,Y(1).signals.values(:,6),'b-','LineWidth',2);
plot(Y(2).time,Y(2).signals.values(:,6),'r--','LineWidth',2);
plot(Y(3).time,Y(3).signals.values(:,6),'k-','LineWidth',2);
plot(Y(4).time,Y(4).signals.values(:,6),'g--','LineWidth',2);
ylim([-1.5,2]);
xlabel('$t$(s)','Interpreter','latex','FontSize',FontSize);
ylabel('$y$(m)','Interpreter','latex','FontSize',FontSize);
h=legend('$y_r$','$y_c$(AC in [22])','$y_c$(PPC in [20])','$y_c$(proposed UC)','$y_c$(proposed PPUC)');
set(gca,'FontName','Helvetica','FontSize',FontSize);
set(h,'Interpreter','latex');

%% ex
figure('color','white');
grid on
hold on
plot(X(1).time,X(1).signals.values(:,4),'m-.','LineWidth',2.5);
plot(X(1).time,X(1).signals.values(:,5),'m-','LineWidth',2);
plot(X(1).time,X(1).signals.values(:,7),'b-','LineWidth',2);
plot(X(2).time,X(2).signals.values(:,7),'r--','LineWidth',2);
plot(X(3).time,X(3).signals.values(:,7),'k-','LineWidth',2);
plot(X(4).time,X(4).signals.values(:,7),'g--','LineWidth',2);
ylim([-0.6,1]);
xlabel('$t$(s)','Interpreter','latex','FontSize',FontSize);
ylabel('$e_x$(m)','Interpreter','latex','FontSize',FontSize);
h=legend('$\rho_x$','$-\rho_x$','$e_x$(AC in [22])','$e_x$(PPC in [20])','$e_x$(proposed UC)','$e_x$(proposed PPUC)');
set(gca,'FontName','Helvetica','FontSize',FontSize);
set(h,'Interpreter','latex');

%% ey
figure('color','white');
grid on
hold on
plot(Y(1).time,Y(1).signals.values(:,4),'m-.','LineWidth',2);
plot(Y(1).time,Y(1).signals.values(:,5),'m-','LineWidth',2);
plot(Y(1).time,Y(1).signals.values(:,7),'b-','LineWidth',2);
plot(Y(2).time,Y(2).signals.values(:,7),'r--','LineWidth',2);
plot(Y(3).time,Y(3).signals.values(:,7),'k-','LineWidth',2);
plot(Y(4).time,Y(4).signals.values(:,7),'g--','LineWidth',2);
ylim([-0.6,1]);
xlabel('$t$(s)','Interpreter','latex','FontSize',FontSize);
ylabel('$e_y$(m)','Interpreter','latex','FontSize',FontSize);
h=legend('$\rho_y$','$-\rho_y$','$e_y$(AC in [22])','$e_y$(PPC in [20])','$e_y$(proposed UC)','$e_y$(proposed PPUC)');
set(gca,'FontName','Helvetica','FontSize',FontSize);
set(h,'Interpreter','latex');

%% \hat{\Delta}_x
figure('color','white');
grid on
hold on
plot(DH(1).time,DH(1).signals.values(:,1),'b-','LineWidth',2);
plot(DH(2).time,DH(2).signals.values(:,1),'r-','LineWidth',2);
plot(DH(3).time,DH(3).signals.values(:,1),'k-','LineWidth',2);
plot(DH(4).time,DH(4).signals.values(:,1),'g--','LineWidth',2);
xlabel('$t$(s)','Interpreter','latex','FontSize',FontSize);
ylabel('$\hat{\Delta}_x$','Interpreter','latex','FontSize',FontSize);
h=legend('$\hat{\Delta}_x$(AC in [22])','$\hat{\Delta}_x$(PPC in [20])','$\hat{\Delta}_x$(proposed UC)','$\hat{\Delta}_x$(proposed PPUC)');
set(gca,'FontName','Helvetica','FontSize',FontSize);
set(h,'Interpreter','latex');

%% \hat{\Delta}_y
figure('color','white');
grid on
hold on
plot(DH(1).time,DH(1).signals.values(:,2),'b-','LineWidth',2);
plot(DH(2).time,DH(2).signals.values(:,2),'r-','LineWidth',2);
plot(DH(3).time,DH(3).signals.values(:,2),'k-','LineWidth',2);
plot(DH(4).time,DH(4).signals.values(:,2),'g--','LineWidth',2);
xlabel('$t$(s)','Interpreter','latex','FontSize',FontSize);
ylabel('$\hat{\Delta}_y$','Interpreter','latex','FontSize',FontSize);
h=legend('$\hat{\Delta}_y$(AC in [22])','$\hat{\Delta}_y$(PPC in [20])','$\hat{\Delta}_y$(proposed UC)','$\hat{\Delta}_y$(proposed PPUC)');
set(gca,'FontName','Helvetica','FontSize',FontSize);
set(h,'Interpreter','latex');
