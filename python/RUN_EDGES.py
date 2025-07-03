from EDGES import EDGES

# Load or simulate your matrices:
X = ...   # ST data
Y = ...   # shared scRNA-seq data
Z = ...   # unique scRNA-seq data
L1 = ...  # normalized graph Laplacian

# Set hyperparameters
lambda1 = 1e-5
lambda2 = 1e1
theta1 = 1e-1
theta2 = 1e-4
tol = 1e-7
d = 20
iterMax = 500

# Run EDGES
W1, W2, H1, H2, p1, p2 = EDGES(X, Y, Z, L1, lambda1, lambda2, theta1, theta2, tol, d, iterMax)