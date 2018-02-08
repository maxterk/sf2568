clear all
close all
load -ascii lecture_note_copy.txt

myMap=colorcube(256);
myMap(1,:)=1;
figure('pos',[10,10,800,800])
pic=reshape(lecture_note_copy,2048,2048);
% pic=histeq(pic,256);
imagesc(pic)
colormap(myMap)
axis square;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
print('mandelbrot','-dpng');

%%
clear all
close all
load -ascii zoom1.txt

myMap=colorcube(256);
myMap(1,:)=1;
figure('pos',[10,10,800,800])
pic=reshape(zoom1,2048,2048);
% pic=histeq(pic,256);
imagesc(pic)
colormap(myMap)
axis square;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
print('mandelbrot_zoom','-dpng');

%%
clear all
close all
load -ascii zoom2.txt

myMap=colorcube(256);
myMap(1,:)=1;
figure('pos',[10,10,800,800])
pic=reshape(zoom2,2048,2048);
% pic=histeq(pic,256);
imagesc(pic)
colormap(myMap)
axis square;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
print('mandelbrot_zoom2','-dpng');

%%
clear all
close all
load -ascii zoom3.txt

myMap=colorcube(256);
myMap(1,:)=1;
figure('pos',[10,10,800,800])
pic=reshape(zoom3,2048,2048);
% pic=histeq(pic,256);
imagesc(pic)
colormap(myMap)
axis square;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
print('mandelbrot_zoom3','-dpng');
