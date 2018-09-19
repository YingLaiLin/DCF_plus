function [ evaluations ] = eval_score(split,ScoreMatrix,use_cdf,use_prc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load('splits_uniform.mat')
evaluations = cell(3,2);
evaluations{1,1} = 'cdf';
evaluations{2,1} = 'recall';
evaluations{3,1} = 'pres';
if use_cdf
    cdf_rates = cdf(full(splits{split}), ScoreMatrix,100);
    evaluations{1,2} = cdf_rates;
end
if use_prc
    rates = recall(full(splits{split}), ScoreMatrix,100) .*100 ;
    pres = precision(full(splits{split}), ScoreMatrix,100) .* 100;
    evaluations{2,2} = rates;
    evaluations{3,2} = pres;   
end

end

