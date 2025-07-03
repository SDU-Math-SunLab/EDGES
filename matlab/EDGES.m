function [W1, W2, H1, H2, p1, p2] = EDGES(X,Y,Z,L1,lambda1,lambda2,theta1,theta2,tol,d,iterMax)
% <INPUT>
%        X: ST data
%        Y: shared scRNA-seq data
%        Z: unique scRNA-seq data
%        L1: normalized ​​graph Laplacian constraint​​ for spatial structure regularization
%        lambda1:spatial hyper-parameter
%        lambda2:sparsity hyper-parameter
%        theta1:decomposition hyper-parameter balancing reconstruction error and constraints for the ST data
%        theta2:decomposition hyper-parameter for the unique scRNA-seq data
%        tol: stopping tolerance
%        d: number of latent factors
%        iterMax: maximum iterations
%        iter: iterations used
% <OUTPUT>
%        p1: denoised expression matrix
%        p2: predicted expression matrix
%        W1: shared gene factors
%        W2: unique scRNA-seq data​​-specific gene factors
%        H1: spatial cell embeddings
%        H2: scRNA-seq data embeddings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    profile('on');
    runtime = java.lang.Runtime.getRuntime();  % for memory monitoring
    startTime = tic;                           % initialize timer
    
    % record initial memory usage
    initialTotal = runtime.totalMemory();
    initialFree = runtime.freeMemory();
    initialUsed = initialTotal - initialFree;
    
    % matrix dimension definitions
    row_W1 = size(X,1);
    row_W2 = size(Z,1);
    col_H1 = size(X,2);
    col_H2 = size(Y,2);
    
    % random initialization with fixed seed for reproducibility 
    rng(42);
    W1 = rand(row_W1,d); 
    W2 = rand(row_W2,d);
    H1 = rand(d,col_H1); 
    H2 = rand(d,col_H2); 
    
    e1d = ones(1,d);
    e1dt = e1d';
   
    % calculate the initial loss
    [delta_init1,~,~,~] = original_loss(W1,W2,H1,H2,X,Y,Z,L1,lambda1,lambda2,theta1,theta2,d);
    delta2 = delta_init1;

    % initialize convergence tracking arrays
    stop = [];
    X1_lossall = [];
    X2_lossall = [];
    X3_lossall = [];
    total_lossall = [];   
    
       %% multiplicative update algorithm
        for iter = 1:iterMax
            disp(iter);
            
            [X1H1t,X2H2t,W1H1H1t,W1H2H2t] = W1_update(W1,H1,H2,X,Y);
            w1 = W1.*((theta1*X1H1t + X2H2t)./(theta1*W1H1H1t + W1H2H2t));
            
            W1 = w1;
                                                         
            [W1tX1,W1tW1H1,eddH1,H1L1] = H1_update(W1,H1,X,L1,d);
            h1 = H1.*((theta1*W1tX1)./(theta1*W1tW1H1 + lambda2*eddH1 + lambda1*H1L1 + 10^(-8)));       
            H1 = h1;
            
            [X3H2t,W2H2H2t] = W2_update(W2,H2,Z);
            w2 = W2.*((X3H2t)./(W2H2H2t));
            
            W2 = w2;            
             
            [W1tW1H2,W1tX2,W2tX3,W2tW2H2,eddH2] = H2_update(W1,W2,H2,Y,Z,d);
            h2 = H2.*((W1tX2 + theta2*W2tX3)./(W1tW1H2 + theta2*W2tW2H2 + lambda2*eddH2));
        
            H2 = h2;  
         
            H1H1t = H1*H1';
            H2H2t = H2*H2';
            
          %% convergence monitoring
            % Calculate reconstruction errors        
            X1_loss = (norm((X - W1*H1),'fro'))^2;
            X2_loss = (norm((Y - W1*H2),'fro'))^2;
            X3_loss = (norm((Z - W2*H2),'fro'))^2;
            H1_loss = trace(H1*L1*H1');
            H_total_loss = e1d*H1H1t*e1dt + e1d*H2H2t*e1dt;
            
            % loss calculation
            total_loss = theta1*X1_loss + X2_loss + theta2*X3_loss+ lambda1*H1_loss + lambda2*H_total_loss; 
            
            % relative loss change calculation
            stop_control1 = abs((delta2 - total_loss)/(delta_init1 - total_loss)); 
            stop_value = abs(stop_control1);

            stop = [stop;stop_value];
            X1_lossall = [X1_lossall;X1_loss];
            X2_lossall = [X2_lossall;X2_loss];
            X3_lossall = [X3_lossall;X3_loss];
            total_lossall = [total_lossall;total_loss];
            
            % check stopping criteria (relative change < tolerance)
                    if stop_value < tol
                        break;
                    end
                        delta2 = total_loss;
        end
        
       %% performance statistics
        totalTime = toc(startTime);  
    
        % calculate memory usage (in MB)
        finalTotal = runtime.totalMemory();
        finalFree = runtime.freeMemory();
        finalUsed = finalTotal - finalFree;
        memUsageMB = (finalUsed - initialUsed) / 1024^2;  
    
       %% display results
        fprintf('\n----- Runtime Statistics -----\n');
        fprintf('Total time elapsed: %.2f seconds\n', totalTime);
        fprintf('Memory usage: %.2f MB\n\n', memUsageMB);
        
       %% output processing
        p1 = W1*H1;
        p2 = W2*H1;
               
       %% save decomposition results
        csvwrite('st_denoise.csv',p1)
        csvwrite('st_inferred.csv',p2)
        
        % Optional: Save decomposition components
        %csvwrite('H1.csv',H1)
        %csvwrite('W1.csv',W1)
        %csvwrite('H2.csv',H2)
        %csvwrite('W2.csv',W2)

end