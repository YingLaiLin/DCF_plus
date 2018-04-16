 
% IMC(1, 200, 100, 100, 10, 0.2);
load ('splitsUniform.mat');


% rate3 = recall(full(splits{1,1}) , ScoreMatrix3 , 300);
% DCF_gen 1.0 model

valid_data = full(splits{1});


% DCF model and DCF_plus model
fileNames = dir('ScoreMatrix*.mat');
fileSize = size(fileNames);
rates = cell(fileSize);
x = 1:100;
hwait = waitbar(0, 'Please wait for data preprocessing');
for iFileIndex = 1:fileSize(1)
    load(fileNames(iFileIndex).name)
    rates{iFileIndex} = cdf(valid_data,ScoreMatrix, 100);    
    str = ['正在运行中', num2str(iFileIndex / fileSize(1) * 100),'%'];
    disp('--------------------------------------')
    disp(fileNames(iFileIndex).name)
    disp(['best rate = ', num2str(rates{iFileIndex}(100))])
    disp(['top 30 avg rate = ', num2str(mean(rates{iFileIndex}(70:100)))]);
    disp(['cdf curves area = ', num2str(trapz(x,rates{iFileIndex}(1:100)))]);
    waitbar(iFileIndex / fileSize(1), hwait, str)
end
close(hwait)

%plot each DCF model with specific name

hold on
for iFileIndex = 1 : fileSize
   
   plot(x, rates{iFileIndex}(1:100), 'linewidth', 2)
end
% plot(x1,y1,'r-',x2,y2,'c-',x3,y3, 'linewidth',2);

lgd = legend(fileNames.name, 'Location', 'best');
lgd.LineWidth = 2;
grid on

% x = 1:100
% plot(x,rate1(1:100),'r-',x,rate2(1:100),'c-','linewidth',2);
% xlabel('Number of genes looked at');
% ylabel('P(hidden gene among genes looked at)');
% grid on

% x = 1:100
% plot(x,pre1(1:100),'r-',x,pre2(1:100),'c-','linewidth',2);
% xlabel('Number of genes looked at');
% ylabel('P(hidden gene among genes looked at)');
% grid on


% clear;
disp('curves complete')
