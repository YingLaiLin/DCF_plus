 
% IMC(1, 200, 100, 100, 10, 0.2);
load ('splitsUniform.mat');


% rate3 = recall(full(splits{1,1}) , ScoreMatrix3 , 300);
% DCF_gen 1.0 model

valid_data = full(splits{1,1});


% DCF model and DCF_gen model, alpha_gen range from 0 to 1, with 0.1 as step
fileNames = dir('ScoreMatrix*.mat');
fileSize = size(fileNames);
precisions = cell(fileSize);
rates = cell(fileSize);
for iFileIndex = 1:fileSize(1)
    load(fileNames(iFileIndex).name)
    rates{iFileIndex} = recall(valid_data , ScoreMatrix, 100) .*100;
    precisions{iFileIndex} = precision(valid_data, ScoreMatrix, 100) .* 100;
    
end

%plot each DCF model with specific name
x = 1:100
hold on
for iFileIndex = 1 : fileSize
   plot(rates{iFileIndex}, precisions{iFileIndex})
end
% plot(x1,y1,'r-',x2,y2,'c-',x3,y3, 'linewidth',2);

legend(fileNames.name)
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
