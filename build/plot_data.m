close all;
load log_file.txt;

figure;
hold;
plot(log_file(:,1), log_file(:,2), 'b');
plot(log_file(:,1), log_file(:,3), 'r');

legend('true', 'estimate');