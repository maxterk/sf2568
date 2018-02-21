close all
clear all
fileID = fopen('output.txt','r');
u=fscanf(fileID,'%f');
u=u';

x=linspace(0,1,numel(u));
h=x(2)-x(1);

% Selected function f
f=sin(5*x(2:end-1));
% Selected function r
r=-exp(x);

ubis=u(3:end)-2*u(2:end-1)+u(1:end-2);
ubis=ubis/h^2;

figure()
plot(x,u);
ru=r.*u;
rhs=ubis+ru(2:end-1);
legend('Solution')
print('Solution','-dpng')
%%
figure();
plot(x(2:end-1),rhs);
hold on
plot(x(2:end-1),f)
fclose('all');

legend('Approximate RHS', 'True, LHS')
print('Error','-dpng')