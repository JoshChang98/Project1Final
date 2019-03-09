function  [mu, Q] = PCA(returns, varargin)

    % Use this function to perform the principal component analysis. After
    % performing the PCA, use the principal components and the eigenvalues
    % to construct a factor model and estimate mu and Q.
    %
    % Note: PCA only requires the asset returns, it does not need the
    % factor returns.

    %----------------------------------------------------------------------

    cov_mat = cov(returns); % creates a covariance matrix of returns
    [V,D] = eig(cov_mat); % V -> matrix of eigenvectors, D -> diag(eigenvalues corresponding to vectors in V)
    [m,n] = size(returns); % gets how big f should be if we don't have K
    
     K=3; 

    
    f = ones(m, K);
    
    
    D=diag(sort(diag(D),'descend')); % make diagonal matrix out of sorted diagonal values of input D
    [c, ind]=sort(diag(D),'descend'); % store the indices of which columns the sorted eigenvalues come from
    V = V(:,ind);
    V
   
    prime = pca(returns)
     diff = flip(V,2) - prime
    V_inv = inv(V);
    for r = 1:m
      for c = 1:K
        f(r, c) =returns(r, :)*V_inv(:, c);
      end
    end
    csvwrite('PCAOutput.csv',V);
    csvwrite('MatlabPCAOutput.csv',prime);
    [mu, Q] = FF(returns, f);

    %----------------------------------------------------------------------

end