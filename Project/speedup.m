%----main program for assignment 1.8-----%

clear all
clc
close all
filename1 = 'output_par_iter_1.txt';
filename2 = 'output_par_iter_2.txt';
filename3 = 'output_par_iter_2.txt';
A1 = importdata(filename1);
A2 = importdata(filename2);
A3 = importdata(filename3);

filename4 = 'output_seq_iter_1.txt';
filename5 = 'output_seq_iter_2.txt';
filename6 = 'output_seq_iter_2.txt';
B1 = importdata(filename4);
B2 = importdata(filename5);
B3 = importdata(filename6);


time_par=(A1+A2+A3)/3;
figure(1)
plot(time_par,'-o')
xlabel('processors')
ylabel('time')

speed_par=time_par(1)./time_par;

figure(2)
plot(speed_par,'-o')
xlabel('processors')
ylabel('speedup')



time_seq=(B1+B2+B3)/3;
figure(3)
plot(time_seq,'-or')
xlabel('processors')
ylabel('time')

speed_seq=time_seq(1)./time_seq;

figure(2)
plot(speed_seq,'-or')
xlabel('processors')
ylabel('speedup')
