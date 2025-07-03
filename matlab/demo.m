% <INPUT>
%        L1: normalized ​​graph Laplacian constraint​​ for spatial structure regularization
%        lambda1:spatial hyper-parameter controlling the weight of the graph Laplacian regularization
%        lambda2:sparsity hyper-parameter enforcing sparse cell representations
%        theta1:decomposition hyper-parameter balancing reconstruction error and constraints for the ST data
%        theta2:decomposition hyper-parameter for the unique scRNA-seq data
%        tol: stopping tolerance
%        d: number of latent factors
%        iterMax: maximum iterations
%        iter: iterations used
%        X1: ST data
%        X2: shared scRNA-seq data
%        X3: unique scRNA-seq data
% <OUTPUT>
%        p1:denoised expression matrix
%        p2: predicted expression matrix
%        W1: shared gene factorss
%        W2: unique scRNA-seq data​​-specific gene factors
%        H1: spatial cell embeddings
%        H2: scRNA-seq data embeddings
%% spatial network construction
 L1 = csvread('laplacian.csv',1,1);

%% hyper-parameter configuration
 lambda1 = 10^(-5);  
 lambda2 = 10^(1); 

 theta1 = 10^(-1);
 theta2 = 10^(-4);

 tol = 10^(-7); 
 d = 20;  
 iterMax = 500;  

%% data loading and preprocessing
% load spatial transcriptomics datasets
% ensure row names of X1 and X2 are in the same order
 X1 = readmatrix('stdata.csv','OutputType', 'string');  
 X2 = readmatrix('shared_scRNA-seq_data.csv','OutputType','string');  
 X3 = readmatrix('unique_scRNAseq_data.csv','OutputType','string');   
 
% align gene names between ST and scRNA-seq data by sorting (must)
 X1 = sortrows(X1, 1);
 X2 = sortrows(X2, 1);
 
% string matrices to numeric
 X1 = str2double(X1);
 X2 = str2double(X2);
 X3 = str2double(X3);
 
% remove the first column (gene names) after alignment
 X1(:,1) = []; 
 X2(:,1) = [];
 X3(:,1) = [];
 
%% EDGES algorithm execution
% Input preprocessed ST data and scRNA-seq data
% The gene order in each fold of X1 must be identical to that in X2 (must)
X11 = csvread('stdata_fold1.csv'); 
X12 = csvread('stdata_fold2.csv');
X13 = csvread('stdata_fold3.csv');

X21 = csvread('shared_scdata_fold1.csv'); 
X22 = csvread('shared_scdata_fold2.csv');
X23 = csvread('shared_scdata_fold3.csv');

X31 = csvread('unique_scdata_fold1.csv'); 
X32 = csvread('unique_scdata_fold2.csv');
X33 = csvread('unique_scdata_fold3.csv');

% Input the validation set for each fold 
% Ensure that the gene order in rows 2001 to the end of X3i is identical to that in X1ipre
X11pre = csvread('stdata_fold1_to_pre.csv'); 
X12pre = csvread('stdata_fold2_to_pre.csv');
X13pre = csvread('stdata_fold3_to_pre.csv');

% prepare data groups for cross-validation
 X_list = {X11, X12, X13}; % Training sets (ST data)
 Y_list = {X21, X22, X23}; % Training sets (shared scRNA-seq data)
 Z_list = {X31, X32, X33}; % Training sets (unique scRNA-seq data)
 
 
% initialize cell array to store results from each fold
 st_inferr_fold = cell(1,3);
 
    for i = 1:3
        X = X_list{i};
        Y = Y_list{i};
        Z = Z_list{i};
        fprintf('Processing the %d th group of data...\n', i);
        
        % execute EDGES algorithm 
         [W1, W2, H1, H2, p1, p2] = EDGES(X, Y, Z, L1, lambda1, lambda2, theta1, theta2, tol, d, iterMax);
   
        % validate row count of p2
         if size(p2, 1) < 2001
             error('p2 in fold %d has less than 2001 rows', i);
         end
    
        % Extract rows from 2001 to end and save
         st_inferr_fold{i} = p2(2001:end, :);  % save inference results
    
    end
    
%% concatenate results from all folds
% ensure the order of st_inferr_fold{i} matches X1ipre
st_infer = vertcat(st_inferr_fold{:}); % predicted data
st_raw = [X11pre;X12pre;X13pre]; % raw stdata
pcc = diag(corr(st_infer',st_raw')); % compute Pearson correlation between inferred and raw data