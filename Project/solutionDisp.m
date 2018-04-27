close all
clear all
fileID = fopen('output.txt','r');
u=fscanf(fileID,'%f');
data=reshape(u,3,numel(u)/3)';

x=data(:,1);
y=data(:,2);
mark=data(:,3)==1;

scatter(x(mark),y(mark),'k','filled')
hold on
scatter(x(~mark),y(~mark),'r','filled')

fclose('all');