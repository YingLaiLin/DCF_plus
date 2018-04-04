% IMCPS(1, 200, 100, 100, 10, 0.2);
load('splitsUniform.mat');
load('genesPhenes.mat');
% load('ScoreMatrix_DCF.mat');

known = full(splits{1});
split = 1;
GenePheneTraining = GenePhene{1};%% 
% oneassoc=full(GenePheneTraining);
numPhenes = size(splits{split},2);  %% 
GenePheneTraining = GenePheneTraining(:,1:numPhenes) - splits{split};%% 
train = full(GenePheneTraining);

% ip = [];
% iq = [];
% id = [];
% %singleton gene
% j = 1;
% for i = 1 : 12331
%     if length(find(known(i,:) == 1)) == 1
%         ip(j,1) = i;
%         j = j + 1;
%     end
% end
% 
% n = 1;
% for m = 1 : 12331
%     if isempty(find(train(m,:) == 1))
%         iq(n,1) = m;
%         n = n + 1;
%     end
% end
% singleton disease
j = 1;
for i = 1 : 3209
    if length(find(known(:,i) == 1)) == 1
        ip(j,1) = i;
        j = j + 1;
    end
end

n = 1;
for m = 1 : 3209
    if isempty(find(train(:,m) == 1))
        iq(n,1) = m;
        n = n + 1;
    end
end

id = intersect(ip,iq);
% save id.mat id;

% known = known(id(:,1),:);
% ScoreMatrix = ScoreMatrix(id(:,1),:);
% rate = recall(known', ScoreMatrix', 100);
% rate = recall(known(id(:,1),:) , ScoreMatrix(id(:,1),:), 100);
load('ScoreMatrix.mat');
rate_DCF = recall(known(:,id(:,1)) , ScoreMatrix(:,id(:,1)), 100);
% precision_IMC = precision(known(:,id(:,1)), ScoreMatrix(:,id(:,1)), 100);
load('ScoreMatrix_0.8.mat');
rate_DCF_8 = recall(known(:,id(:,1)) , ScoreMatrix(:,id(:,1)), 100);
% DCF+0.8
x = 1:100
figure();
plot(x,rate_DCF(1:100),'-',x, rate_DCF_8(1:100),'--');

xlabel('Number of genes looked at');
ylabel('P(hidden gene among genes looked at)');
grid on

% save rate_DCF1.mat rate_DCF1;

% save singletonrate_IMC1.mat rate_IMC1;
% save singletonpre_IMC.mat precision_IMC;


% x = rate_IMC1 .* 100;
% y = precision_IMC .* 100;
% figure(2);
% plot(x,y,'y-','linewidth',2);
% % plot(x1,y1,'c--','linewidth',2);
% xlabel('Recall(%)');
% ylabel('Precision(%)');
% grid on




