import numpy as np
from scipy.linalg import norm

def EDGES(X, Y, Z, L1, lambda1, lambda2, theta1, theta2, tol, d, iterMax):
    """
    EDGES: Joint decomposition model for spatial transcriptomics and scRNA-seq integration.

    Parameters:
        X      : ST data matrix (shared genes × spatial cells)
        Y      : Shared scRNA-seq data matrix (shared genes × single cells)
        Z      : Unique scRNA-seq data matrix (unique genes × single cells)
        L1     : Normalized graph Laplacian (spatial adjacency constraint)
        lambda1: Graph Laplacian regularization weight
        lambda2: Sparsity regularization weight
        theta1 : Weight for ST reconstruction loss
        theta2 : Weight for unique scRNA-seq loss
        tol    : Tolerance for convergence
        d      : Number of latent factors
        iterMax: Maximum number of iterations

    Returns:
        W1: Shared gene factor matrix
        W2: Unique gene factor matrix (scRNA-seq-specific)
        H1: ST cell embeddings
        H2: scRNA-seq cell embeddings
        p1: Denoised ST data matrix
        p2: Predicted expression matrix for ST cells
    """
    # Initialize matrix shapes
    row_W1 = X.shape[0]
    row_W2 = Z.shape[0]
    col_H1 = X.shape[1]
    col_H2 = Y.shape[1]
    
    # Random initialization with fixed seed 
    rng = np.random.RandomState(seed=42)
    W1 = rng.rand(row_W1, d)
    W2 = rng.rand(row_W2, d)
    H1 = rng.rand(d, col_H1)
    H2 = rng.rand(d, col_H2)

    e1d = np.ones((1, d))

    # Compute initial loss before starting iterations
    loss_X1 = norm(X - W1 @ H1, 'fro')**2
    loss_X2 = norm(Y - W1 @ H2, 'fro')**2
    loss_X3 = norm(Z - W2 @ H2, 'fro')**2
    loss_H1 = np.trace(H1 @ L1 @ H1.T)
    loss_H = e1d @ (H1 @ H1.T) @ e1d.T + e1d @ (H2 @ H2.T) @ e1d.T

    delta_init1 = theta1 * loss_X1 + loss_X2 + theta2 * loss_X3 + lambda1 * loss_H1 + lambda2 * loss_H
    delta2 = delta_init1

    # Main optimization loop
    for iter in range(iterMax):
        # Update W1
        X1H1t = X @ H1.T
        X2H2t = Y @ H2.T
        W1H1H1t = W1 @ (H1 @ H1.T)
        W1H2H2t = W1 @ (H2 @ H2.T)
        W1 *= (theta1 * X1H1t + X2H2t) / (theta1 * W1H1H1t + W1H2H2t + 1e-8)

        # Update H1
        W1tX1 = W1.T @ X
        W1tW1H1 = W1.T @ W1 @ H1
        eddH1 = (e1d.T @ e1d) @ H1
        H1L1 = H1 @ L1
        H1 *= (theta1 * W1tX1) / (theta1 * W1tW1H1 + lambda2 * eddH1 + lambda1 * H1L1 + 1e-8)

        # Update W2
        X3H2t = Z @ H2.T
        W2H2H2t = W2 @ (H2 @ H2.T)
        W2 *= X3H2t / (W2H2H2t + 1e-8)

        # Update H2
        W1tX2 = W1.T @ Y
        W2tX3 = W2.T @ Z
        W1tW1H2 = W1.T @ W1 @ H2
        W2tW2H2 = W2.T @ W2 @ H2
        eddH2 = (e1d.T @ e1d) @ H2
        H2 *= (W1tX2 + theta2 * W2tX3) / (W1tW1H2 + theta2 * W2tW2H2 + lambda2 * eddH2 + 1e-8)

        # Compute reconstruction loss
        loss_X1 = norm(X - W1 @ H1, 'fro')**2
        loss_X2 = norm(Y - W1 @ H2, 'fro')**2
        loss_X3 = norm(Z - W2 @ H2, 'fro')**2
        loss_H1 = np.trace(H1 @ L1 @ H1.T)
        loss_H = e1d @ (H1 @ H1.T) @ e1d.T + e1d @ (H2 @ H2.T) @ e1d.T

        total_loss = (theta1 * loss_X1 + loss_X2 + theta2 * loss_X3 +
                      lambda1 * loss_H1 + lambda2 * loss_H)

        # Convergence check using relative change in loss
        stop_value = abs((delta2 - total_loss) / (delta_init1 - total_loss + 1e-8))
        if stop_value < tol:
            print(f"Converged at iteration {iter+1}")
            break

        delta2 = total_loss

    # Compute final outputs
    p1 = W1 @ H1  # Denoised ST matrix
    p2 = W2 @ H1  # Predicted expression matrix from unique scRNA-seq

    return W1, W2, H1, H2, p1, p2
