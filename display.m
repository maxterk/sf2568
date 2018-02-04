clear all,close all
load -ascii color.txt

 
 color_array=zeros(2048,2048);
 for n=1:2048
     for m=1:2048
         color_array(n,m)=color(m +(n-1)*2048);
     end
 end
 

imagesc(color_array')
colormap(colorcube)
