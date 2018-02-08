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
load -ascii color.txt

myMap=colorcube(256);
myMap(1,:)=1;
figure('pos',[10,10,800,800])
pic=reshape(color,2048,2048);
% pic=histeq(pic,256);
imagesc(pic)
colormap(myMap)
axis square;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
print('mandelbrot_zoom','-dpng');
