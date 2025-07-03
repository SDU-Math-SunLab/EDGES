function [delta_init1,X1_loss_raw,X2_loss_raw,X3_loss_raw] = original_loss(W1,W2,H1,H2,X,Y,Z,L1,lambda1,lambda2,theta1,theta2,d)
% <INPUT>
%        W1: initial shared gene factors
%        W2: initial ​​unique scRNA-seq data​​-specific gene factors
%        H1: initial spatial cell embeddings
%        H2: initial ​​scRNA-seq data embeddings
%        X: ST data
%        Y: shared scRNA-seq data
%        Z: unique scRNA-seq data
%        L1: normalized ​​graph Laplacian constraint​​ for spatial structure regularization
%        lambda1:spatial hyper-parameter controlling the weight of the graph Laplacian regularization
%        lambda2:sparsity hyper-parameter enforcing sparse cell representations
%        theta1:decomposition hyper-parameter balancing reconstruction error and constraints for the ST data
%        theta2:decomposition hyper-parameter for the unique scRNA-seq data
%        tol: stopping tolerance
%        d: number of latent factors
%        iterMax: maximum iterations
%        iter: #iterations used
% <OUTPUT>
%        delta_init1: initial ​​composite loss
%        X1_loss_raw: initial ​​ST reconstruction loss​​
%        X2_loss_raw: initial ​​shared scRNA-seq data reconstruction loss​​
%        X3_loss_raw: initial ​​unique scRNA-seq data reconstruction loss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H1t = H1';
    H2t = H2';

    H2H2t = H2*H2t;
    H1H1t = H1*H1t;
    
    e1d = ones(1,d);
    e1dt = e1d';
    eH1H1tet = e1d*H1H1t*e1dt;
    eH2H2tet = e1d*H2H2t*e1dt;
    
    X1_loss_raw = (norm((X - W1*H1),'fro'))^2;
    X2_loss_raw = (norm((Y - W1*H2),'fro'))^2;
    X3_loss_raw = (norm((Z - W2*H2),'fro'))^2;
    H1_loss_raw = trace(H1*L1*H1');
    H_totalloss_raw = eH1H1tet + eH2H2tet;
    
    delta_init1 = theta1*X1_loss_raw + X2_loss_raw + theta2*X3_loss_raw + lambda1*H1_loss_raw + lambda2*H_totalloss_raw; % # all loss

end



                