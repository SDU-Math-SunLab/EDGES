function [W1tX1, W1tW1H1, eddH1, H1L1] = H1_update(W1, H1, X, L1, d)
% <INPUT>
%        W1: shared gene factors
%        H1: spatial cell embeddings
%        X: ST data
%        L1: normalized ​​graph Laplacian constraint​​ for spatial structure regularization
%        d: number of latent factors
% <OUTPUT>
%        W1tX1
%        W1tW1H1
%        ekkH1
%        H1L1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    W1t = W1';
    H1L1 = H1*L1;
    W1tX1 = W1t*X;
    W1tW1H1 = W1t*W1*H1;
    edd = ones(d,d);
    eddH1 = edd*H1;   
    
end