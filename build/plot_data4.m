close all;
load log_file.txt;


hold;
plot(log_file(:,1), log_file(:,2), 'b');
plot(log_file(:,1), log_file(:,3), '--r');
xlabel('time');
ylabel('f(x, )');
title('position');


legend('true', 'estimate');