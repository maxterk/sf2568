clear all
close all
load -ascii color.txt

myMap=colorcube(256);
myMap(1,:)=1;
figure('pos',[10,10,800,800])
pic=reshape(color,2048,2048);
% pic=histeq(pic,256);
imagesc(pic)
colormap(myMap)
axis square;
print('mandelbrot','-dpng');
