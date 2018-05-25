
function ScoreMatrix = IMC(split, k, featureRank, networkRank, loss, lambda, emb_dim)
%%%% INPUT PARAMS:
%%% split       : One of the k-fold splits to be hidden.  %% 留一验证法, 记住用于训练的集合
%%% k           : Rank of the parameter matrix Z. %% 作为参数的矩阵的秩
%%% featureRank : Number of compressed dimensions (microarray gene expression data). %% 希望得到的数据降维之后的维度
%%% networkRank : Number of latent features from similarity networks. %% 相似性网络的维度
%%% loss	    : loss function to be used (10 - for squared loss) %% 损失函数的参数
%%% lambda      : Regularization parameter l in the objective. %% ?
%%% alpha       : the parameter r for PU formulation  %% ?

%%%% OUTPUT:
%%% ScoreMatrix : Matrix of size # genes x # diseases; Each column gives the scores estimated by IMC for the corresponding OMIM disease.	 %% 评分矩阵

%%% Data files are provided in the same directory.

% load data from file
genesPhenes = load ('genesPhenes.mat');  %% 基因表型
geneFeatures = load ('geneFeatures.mat');	%% 基因特征
clinicalFeatures = load ('clinicalFeaturesTFidf.mat'); %% 临床特征

% specify filenames to save .mat files
REDUCTION_MICROARRAY_FILENAME = 'reduction_microarray';
REDUCTION_PHENOTYPES_FILENAME = 'reduction_phenotypes';
POSTFIX_OF_DATA = '.mat';

% control whether to execute embedding
can_phenetype_embedding = true;
can_imfTrain = false;

% 取第一个物种基因集作为训练集
genePheneTraining = genesPhenes.GenePhene{1};
if (split > 0) %% else use all data to m ake predictions.
    splits = load('splitsUniform.mat');
    numPhenes = size(splits.splits{split},2); %% 选择其中一个 cell 作为测试集, 并返回 splits 列数
    % 由于矩阵的值都是1, 故两个矩阵相减就是把一个矩阵作为训练集, 另外一个作为验证集
    genePheneTraining = genePheneTraining(:,1:numPhenes) - splits.splits{split}; 
    
end
features = [];  %% 行特征 -- gene
colFeatures = []; %% 列特征 -- disease

% 分别对 Gene 和 Disease 进行降维
if (featureRank == Inf)	%%如果需要的维度很大, 那么不需要降维, 直接返回默认矩阵
    features = geneFeatures.GeneFeatures;  
%     colFeatures = clinicalFeatures.colFeatures;	
    colFeatures = clinicalFeatures.F;	
elseif (featureRank > 0)
    disp ('Reducing feature dimensionality..');
    % Reducing dimensionality of microarray
    [U,S] = svd(geneFeatures.GeneFeatures' * geneFeatures.GeneFeatures);
    features = geneFeatures.GeneFeatures * U(:, 1 : featureRank);
    % Reducing dimensionality of OMIM pages
    [U,S] = svds(sparse(clinicalFeatures.F(1:numPhenes,:)),featureRank);
    colFeatures = U;

end

if (networkRank > 0)
    disp ('Reducing gene network dimensionality..');
    % Reducing dimensionality of gene-gene interaction network
    %     mysdae(full(genesPhenes.GeneGene_Hs), networkRank, 1);
    %     sdae_data = load('sdae999.mat');
    %     features = [features sdae_data.H(1:gene_phenes.numGenes,:)];
    %     clear sdae_data.H
    network_filename = 'HumanNet_mapped.mat';
    network_features = load(network_filename);
    features = [features network_features.features(1:genesPhenes.numGenes,:) ];
    % Reducing dimensionality of orthologous phenotypes, 循环处理 GenePhene
    % 中的每个 Cell
    %GP_sp = [];
    %for sp=2:numel(genesPhenes.GenePhene)
    %    GP = genesPhenes.GenePhene{sp};
    %    for i=1:size(GP,2)
    %        if(sum(GP(:,i)) > 0)
    %            GP(:,i) = GP(:,i)/norm(GP(:,i));
    %        end
    %    end
    %    GP_sp = [GP_sp GP];
    %end
    % 同源基因
%     phenotypes_filename = [REDUCTION_PHENOTYPES_FILENAME,POSTFIX_OF_DATA];
%     if can_phenetype_embedding || ~exist(phenotypes_filename,'file')
%         mysdae(full(GP_sp), networkRank, 1, REDUCTION_PHENOTYPES_FILENAME);
%     end
%     sdae_data = load(phenotypes_filename);
%     features = [features sdae_data.H];
%     clear sdae_data.H
%     
    % Reducing dimensionality of  disease similarities network
    uFilename = 'U.mat';
    if ~exist(uFilename,'file')
        [U,S] = svds(genesPhenes.PhenotypeSimilaritiesLog(1:numPhenes,1:numPhenes), networkRank);
        save uFilename U
        clear S
    else
        load(uFilename);
    end
    colFeatures = [colFeatures U];
    clear U
    
end

for i=1:size(features,2)
    if (norm(features(:,i)) > 0)
        features(:,i) = features(:,i)/norm(features(:,i)); %% 返回矢量或者矩阵的2-范数
    end
end
for i=1:size(colFeatures,2)
    if (norm(colFeatures(:,i))>0)
        colFeatures(:,i) = colFeatures(:,i)/norm(colFeatures(:,i));
    end
end

% If no features are supplied (i.e. networkRank & featureRank are 0), use default matrix completion.
if (isempty(features))
    features = eye(genesPhenes.numGenes);
end
if (isempty(colFeatures))
    colFeatures = eye(numPhenes);
end
fprintf('Using %d row, %d col features.\n', size(features,2), size(colFeatures,2));

matrix_filename = 'P.mat';
if can_imfTrain || ~exist(matrix_filename,'file')
    [W H time] = train_mf(sparse(double(genePheneTraining)), sparse(double(features)), sparse(double(colFeatures)), [' -l ' num2str(lambda) ' -k ' num2str(k) ' -t 10' ' -s ' num2str(loss)]); 
    save 'P.mat' W H ;
else
    load(matrix_filename);
end

ScoreMatrix = features * W *H' * colFeatures';
scoreMatrixFilename = sprintf('IMC_ScoreMatrix_lambda_%.2f_dim_a%d.mat',lambda,emb_dim);
save(scoreMatrixFilename,'ScoreMatrix');
send_mail_upon_finished('HN2vec ScoreMatrix got',strcat(scoreMatrixFilename,' got') , '18850544602@163.com');
end


