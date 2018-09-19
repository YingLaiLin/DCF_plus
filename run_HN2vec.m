%filenames = dir('6on10*.mat');
%s = size(filenames);
%disp(s);
%for findex = 1:s(1)
%    HN2vec(1,200,100,100,0.01,0.01,filenames(findex).name);
%end
HN2vec(1,200,100,100,0.01,0.01,'GeneGeneHsp2q1.mat');
exit;
