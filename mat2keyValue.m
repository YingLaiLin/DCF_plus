load('genesPhenes.mat')
fout = fopen('ggi_weighted.txt', 'w');
[row, col, v] = find(GeneGene_Hs);
nData = size(row);
genes = 12331;
total_values = nData(1)
fprintf(fout, "%d %d\n", genes, total_values);
formatSpec = "%d %d %f\n";

for iData = 1:total_values
    if mod(iData,10000) == 0
        disp('10000 data completed')
    end
    fprintf(fout, formatSpec, row(iData), col(iData), v(iData));
    
end
back = fclose(fout);
disp('FINISH')
