function [X3H2t,W2H2H2t] = W2_update(W2,H2,Z)
% <INPUT>
%        W2: unique scRNA-seq gene factors
%        H2: scRNA-seq cell embeddings
%        Z: unique scRNA-seq data
% <OUTPUT>
%        X3H2t
%        W2H2H2t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    H2t = H2';   
    X3H2t = Z*H2t;
    H2H2t = H2*H2t;
    W2H2H2t = W2*H2H2t;
  
end