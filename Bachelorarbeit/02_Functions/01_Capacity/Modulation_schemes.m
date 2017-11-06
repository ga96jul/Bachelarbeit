x = [1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1 3 3 3 3 3 3 3 3 -3 -3 -3 -3 -3 -3 -3 -3 5 5 5 5 5 5 5 5 -5 -5 -5 -5 -5 -5 -5 -5 7 7 7 7 7 7 7 7 -7 -7 -7 -7 -7 -7 -7 -7];
y = [7 5 3 1 -1 -3 -5 -7 7 5 3 1 -1 -3 -5 -7 7 5 3 1 -1 -3 -5 -7 7 5 3 1 -1 -3 -5 -7 7 5 3 1 -1 -3 -5 -7 7 5 3 1 -1 -3 -5 -7 7 5 3 1 -1 -3 -5 -7 7 5 3 1 -1 -3 -5 -7];
figure;
plot(x,y,'*');
hold on;
grid on;
title('64-QAM modulation');
xlabel('I');
ylabel('Q');
xlim([-8, 8]);
ylim([-8, 8]);