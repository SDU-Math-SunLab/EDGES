function [W1tW1H2,W1tX2,W2tX3,W2tW2H2,eddH2] = H2_update(W1,W2,H2,Y,Z,d)
% <INPUT>
%        W1: shared gene factors
%        W2: unique scRNA-seq gene factors
%        H2: scRNA-seq cell embeddings
%        Y: shared scRNA-seq data
%        Z: unique scRNA-seq data
%        d: number of latent factors
% <OUTPUT>
%        W1tW1H2
%        W1tX2
%        W2tX3
%        W2tW2H2
%        ekkH2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    W1t = W1';
    W2t = W2'; 
    W1tX2 = W1t*Y;
    W2tX3 = W2t*Z; 
    W1tW1H2 = W1t*W1*H2;
    W2tW2H2 = W2t*W2*H2; 
    edd = ones(d,d);
    eddH2 = edd*H2;
    
end