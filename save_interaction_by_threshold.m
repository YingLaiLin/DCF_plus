function  save_interaction_by_threshold(threshold, is_used_for_weighted_graph)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load('genesPhenes.mat')
% set filename for different situation
if is_used_for_weighted_graph
    ggi_filename =  sprintf('ggi_%.1f_weighted.txt',threshold);
else
    ggi_filename =  sprintf('ggi_%.1f_unweighted.txt',threshold);
end
fout = fopen(ggi_filename, 'w');
% get triple tuple 
GeneGene_Hs(GeneGene_Hs < threshold) = 0;
[row, col, v] = find(GeneGene_Hs);
nData = size(row);
genes = 12331; % fixed size of gene
total_values = nData(1); % assume dimension is (N,1)
fprintf(fout, "%d %d\n", genes, total_values);
if is_used_for_weighted_graph
    formatSpec = "%d %d %f\n";
    for iData = 1:total_values
        if mod(iData,10000) == 0
            disp('10000 data completed')
        end
        fprintf(fout, formatSpec, row(iData), col(iData), v(iData));    
    end
else
    formatSpec = "%d %d\n";
    for iData = 1:total_values
        if mod(iData,10000) == 0
            disp('10000 data completed')
        end
        fprintf(fout, formatSpec, row(iData), col(iData));
    end
end
fclose(fout);
disp('ggi extraction finish')



