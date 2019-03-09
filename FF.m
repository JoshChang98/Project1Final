function  [mu, Q] = FF(returns, factRet)    
% Use this function to perform a linear regression model for CAPM.
% Number of observations;
N = size(returns, 1);
%Number of columns in returns
n = size(returns,2);
%Number of columns in factor table
nFactCols = size(factRet,2);
% Calculate the factor expected excess return from historical data using
% the geometric mean
expExFactRet = geomean(factRet + 1) - 1;
% Calculate the factor variance
sigmaF = var(factRet);    
%Creating the matrix X of ones and the factor
 X = [ones(N,1) factRet];
%Use the closed-form (CF) solution to find the collection of alphas 
% and betas for all assets
temp = inv(transpose(X)* X)*transpose(X)*returns;
alpha  = temp(1,:);
%Creating betaCF
nTempRows = size(temp,1);
nTempCols = size(temp,2);
temp_betaCF = zeros(nTempRows-1,nTempCols);
for i = 2: nTempRows
    temp_betaCF(i-1,:) = temp(i,:);
end
betaCF = transpose(temp_betaCF);   

%Calculating matrix of sum of intercepts (alphas) and factors multiplied by
%factor loading
temp2 = X*temp;

%Calculating residual matrix
resid = returns - temp2;
denom = N -nFactCols-1;
%Calculating diagonal matrix of residuals
D = zeros(n);
    for i = 1:n
       D(i,i) = sumsqr(resid(:,i))/ denom; 
    end

    
% %Calculate matrix of factor loading multiplied by expected excess factor
Beta_x_Fact = zeros(nTempRows-1,nTempCols);
for i = 1:nFactCols
   Beta_x_Fact(i,:) = expExFactRet(i) * temp(i+1,:); 
end
sum_Beta_x_Fact = sum(Beta_x_Fact,1);
mu = transpose(alpha) + transpose(sum_Beta_x_Fact)
%temp2
temp3 = geomean(temp2+1,1)-1;
mu_2 = transpose(temp3)
%Calculating the factor covariance matrix
F = cov(factRet);
%Calculating the asset covariance matrix
Q = betaCF *F*transpose(betaCF) + D;
    % mu =          % n x 1 vector of asset exp. returns
    % Q  =          % n x n asset covariance matrix
    %----------------------------------------------------------------------
    
end