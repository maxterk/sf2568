close all
clear all
fileID = fopen('output.txt','r');
u=fscanf(fileID,'%f');
u=u';

x=linspace(0,1,numel(u));
h=x(2)-x(1);
x2=x.^2;

% Selected function f
f=sin(5*x2(2:end-1));
% Selected function r
r=x2-x-1;

ubis=u(3:end)-2*u(2:end-1)+u(1:end-2);
ubis=ubis/h^2;

figure()
plot(x,u);
ru=r.*u;
rhs=ubis+ru(2:end-1);
xlabel('x')
ylabel('u')
legend('Solution')
print('Solution','-dpng')
fclose('all');
%%
figure();
plot(x(2:end-1),rhs-f);
xlabel('x')
ylabel('e(x)')

legend('Error')
print('Error','-dpng')