close all;
load log_file.txt;

subplot(2,1,1);
hold;
plot(log_file(:,1), log_file(:,2), 'b');
plot(log_file(:,1), log_file(:,4), 'r');
xlabel('position');

subplot(2,1,2);
hold
plot(log_file(:,1), log_file(:,3), 'b');
plot(log_file(:,1), log_file(:,5), 'r');
xlabel('velocity');

legend('true', 'estimate');