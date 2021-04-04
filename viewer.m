function viewer
clc
close all
A=load('spectrum.txt');
w=70.0e9;



f=figure(1)

hold all
plot(A(:,1)/w,A(:,2)*w,'linewidth',1.5)
plot(A(:,1)/w,A(:,3)*w,'linewidth',1.5)
plot(A(:,1)/w,A(:,4)*w,'linewidth',1.5)
plot(A(:,1)/w,A(:,5)*w,'linewidth',1.5)

plot(A(:,1)/w,A(:,6)*w,'linewidth',1.5)
plot(A(:,1)/w,A(:,7)*w,'linewidth',1.5)
plot(A(:,1)/w,A(:,8)*w,'linewidth',1.5)
plot(A(:,1)/w,A(:,9)*w,'linewidth',1.5)


NumTicks = 5;
L = get(gca,'YLim');
% set(gca,'YTick',linspace(0,0.75e-3,NumTicks))
set(gca,'Fontsize', 14)
xlabel('\epsilon_-/\omega','Interpreter','tex','fontsize',15)
ylabel('\omega dP/d\epsilon_-','Interpreter','tex','fontsize',15)


%lgd = legend('$\uparrow\uparrow\mathbf{e}_\parallel$','$\uparrow\downarrow\mathbf{e}_\parallel$','$\uparrow\downarrow\mathbf{e}_\perp$');
%set(lgd,'Fontsize',13)

f1=figure(2)
set(f1,'DefaultTextInterpreter','Latex')
set(f1,'DefaultAxesTickLabelInterpreter','Latex')
#set(f1,'DefaultlegendInterpreter','Latex')
hold all
plot(A(:,1)/w,w*(A(:,2)+A(:,3)+A(:,4)+A(:,5)),'linewidth',1.5)
plot(A(:,1)/w,w*(A(:,6)+A(:,7)+A(:,8)+A(:,9)),'linewidth',1.5)

NumTicks = 5;
L = get(gca,'YLim');
% set(gca,'YTick',linspace(0,0.75e-3,NumTicks))
set(gca,'Fontsize', 14)
xlabel('\epsilon_-/\omega','Interpreter','tex','fontsize',15)
ylabel('\omega dP/d\epsilon_-','Interpreter','tex','fontsize',15)

