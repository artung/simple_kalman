close all;
load log_file.txt;


hold;
plot(log_file(:,1), log_file(:,2), 'b');
plot(log_file(:,5), log_file(:,6), 'r');
xlabel('x1');
ylabel('x2');
title('position');


legend('true', 'estimate');