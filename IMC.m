% Make predictions for OMIM phenotypes using IMC: Natarajan, N., Dhillon, I. S. Inductive Matrix Completion for Predicting Gene-Disease Associations. To appear in Bioinformatics, 2014. 
% 
% Programmer: Nagarajan Natarajan.
% Last update: Mar 21, 2014.
%

function ScoreMatrix = IMC(split, k, featureRank, networkRank, loss, lambda)
    %%%% INPUT PARAMS:
    %%% split       : One of the k-fold splits to be hidden.
    %%% k           : Rank of the parameter matrix Z.
    %%% featureRank : Number of PCA dimensions (microarray gene expression data).
    %%% networkRank : Number of latent features from similarity networks.
    %%% loss	    : loss function to be used (10 - for squared loss) 
    %%% lambda      : Regularization parameter in the objective.

    %%%% OUTPUT:
    %%% ScoreMatrix : Matrix of size # genes x # diseases; Each column gives the scores estimated by IMC for the corresponding OMIM disease.
    
	%%% Data files are provided in the same directory.
    load ('genesPhenes.mat');
	load ('geneFeatures.mat');
	load ('clinicalFeaturesTFidf.mat');
	GenePheneTraining = GenePhene{1};
	if (split > 0) %% else use all data to make predictions.
	    load ('splitsUniform.mat');
	    numPhenes = size(splits{split},2);
	    GenePheneTraining = GenePheneTraining(:,1:numPhenes) - splits{split};
	end
	Features = [];
	ColFeatures = [];
	if (featureRank == Inf)
		Features = GeneFeatures;
		ColFeatures = F;
	elseif (featureRank > 0)
		disp ('Reducing feature dimensionality..');
		[U,~] = svd(GeneFeatures'*GeneFeatures);
		Features = GeneFeatures * U(:,1:featureRank);
		[U,~] = svds(sparse(F(1:numPhenes,:)), featureRank);
		ColFeatures = U;
	end
	if (networkRank > 0)
		disp ('Reducing gene network dimensionality..');
		[U,S] = svds(GeneGene_Hs, networkRank);
		Features = [Features U(1:numGenes,:)];
		clear U S
		GP_sp = [];
    		for sp=2:numel(GenePhene)
        		GP = GenePhene{sp};
        		for i=1:size(GP,2)
            		if(sum(GP(:,i)) > 0)
                		GP(:,i) = GP(:,i)/norm(GP(:,i));
            		end
        		end
			GP_sp = [GP_sp GP];
		end
		[U,S] = svds(GP_sp, networkRank);
		Features = [Features U];
		clear U S
		[U,~] = svds(PhenotypeSimilaritiesLog(1:numPhenes,1:numPhenes), networkRank);
		ColFeatures = [ColFeatures U];
		clear U S
	end
	for i=1:size(Features,2)
		if (norm(Features(:,i)) > 0)
			Features(:,i) = Features(:,i)/norm(Features(:,i));
		end
	end
	for i=1:size(ColFeatures,2)
		if (norm(ColFeatures(:,i))>0)
			ColFeatures(:,i) = ColFeatures(:,i)/norm(ColFeatures(:,i));
		end
	end
        %% If no features are supplied (i.e. networkRank & featureRank are 0), use default matrix completion.
	if (isempty(Features))
		Features = eye(numGenes);
	end
	if (isempty(ColFeatures))
		ColFeatures = eye(numPhenes);
	end
	fprintf('Using %d row, %d col features.\n', size(Features,2), size(ColFeatures,2));
	%% Invoke train_mf function from the IMC library. 
	[W, H, ~] = train_mf(sparse(GenePheneTraining), sparse(Features), sparse(ColFeatures), [' -l ' num2str(lambda) ' -k ' num2str(k) ' -t 10' ' -s ' num2str(loss)]); 
	ScoreMatrix = Features * W'*H * ColFeatures';
    scoreMatrixFilename = sprintf('IMC_ScoreMatrix_lambda_%.2f.mat',lambda);
    save(scoreMatrixFilename,'ScoreMatrix');
%    load('splitsUniform.mat');
%    cdf_rates = cdf(full(splits{1}),ScoreMatrix,100);
%    rates = recall(full(splits{1}),ScoreMatrix, 100) .* 100;
%    pres = precision(full(splits{1}),ScoreMatrix,100) .* 100;
%    send_mail_upon_finished('IMC ScoreMatrix got',sprintf('cdf best score=%.4f AUPRC=%.4f', cdf_rates(100), trapz(rates(1:100),pres(1:100))), '18850544602@163.com');
%    save(name, 'cdf_rates','rates','pres');
end


