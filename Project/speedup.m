clear all
clc
close all
filename1 = 'output_20000_synch_iter_2.txt';
filename2 = 'output_20000_synch_iter_3.txt';
filename3 = 'output_20000_synch_iter_1.txt';
A1 = importdata(filename1);
A2 = importdata(filename2);
A3 = importdata(filename3);

filename4 = 'output_20000_non_synch_iter_1.txt';
filename5 = 'output_20000_non_synch_iter_2.txt';
filename6 = 'output_20000_non_synch_iter_3.txt';
B1 = importdata(filename4);
B2 = importdata(filename5);
B3 = importdata(filename6);

time_par=(A1+A2+A3)/3;
time_seq=(B1+B2+B3)/3;

%%
np=(1:numel(time_seq))';
cumSeq=np.*time_seq;
figure(1)
plot(cumSeq,'-o')
fastestSeq=min(cumSeq);
xlabel('Segments')
ylabel('Sequential time')
legend('No synchronization')
print('sequentialDivConq','-dpng')

%%
figure(2)   
plot(time_par,'-o')
hold on
plot(time_seq,'-or')
xlabel('Processors')
ylabel('Iteration time')
legend('Synchronization','No synchronization')
print('parallelDivConq','-dpng')

%%
speed_par=fastestSeq./time_par;
speed_seq=fastestSeq./time_seq;
figure(3)
plot(speed_par,'-o')
hold on
plot(speed_seq,'-or')
plot(1:numel(time_seq))
xlabel('Processors')
ylabel('Speedup')
legend('Synchronization','No synchronization','Linear')
print('speedUp','-dpng')

%%
speed_par=time_par(1)./time_par;
speed_seq=time_seq(1)./time_seq;
figure(4)
plot(speed_par,'-o')
hold on
plot(speed_seq,'-or')
plot(1:numel(time_seq))
xlabel('Processors')
ylabel('Speedup')
legend('Synchronization','No synchronization','Linear')
print('notSpeedUp','-dpng')

%%
figure(5)
plot(100*time_par./time_seq-100)
xlabel('Processors')
ylabel('Communication overhead time [%]')
print('Overhead','-dpng')

%%
speed_par=fastestSeq./time_par;
figure(6)
plot(speed_par./np,'-o')
xlabel('Processors')
ylabel('Efficiency')
legend('Synchronization')
print('Efficiency','-dpng')