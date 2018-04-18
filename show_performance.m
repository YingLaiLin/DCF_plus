load('splitsUniform.mat')
load('ScoreMatrix.mat')
topR = 100;
% rate = recall(full(splits{1,1}), ScoreMatrix, topR);
precision = precision(full(splits{1,1}), ScoreMatrix, topR);
x = 1:100
% plot(x, rate(1:100), '--')
 plot(x, precision(1:100), '--')
xlabel('Number of genes looked at')
ylabel('P(hidden gene amoing genes looked at)')
grid on
