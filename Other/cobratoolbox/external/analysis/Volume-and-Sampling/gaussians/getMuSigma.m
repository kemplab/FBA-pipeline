function [ mu, sigma ] = getMuSigma(model, rP, centroid)
%GETMUSIGMA This function takes in a metabolic network model, and gives a
%mean and covariance for a Gaussian distribution over the
%polytope. The mean and covariance are probably not the best choice, but
%this could be viewed as a first guess.
%Input:

% model: the metabolic model
% rP: the rounded polytope corresponding to model
% centroid: an estimate for the centroid of the polytope, if already known

if nargin < 2
    rP = preprocess(chrrParseModel(model)); 
end

if nargin < 3
   centroid = getCentroid(rP);
end

mu = centroid;

%we'll choose the covariance matrix to be a diagonal matrix
%I am not sure how to infer covariances between different reactions.
numRxns = length(model.rxns);

%only infer covariances of exchange reactions
[model] = findSExRxnInd(model);
model_exch = ~model.SIntRxnBool;
sigma_diag = zeros(numRxns,1);
%build the diagonal of sigma, then convert it to a sparse matrix
for i=1:numRxns
   if model_exch(i)==1
       sigma_diag(i) = ((model.ub(i)-model.lb(i))/400)^2;
   else
       %if not an exchange reaction, set it to something big
       %note: doing uniform sampling is equivalent to Gaussian
       %sampling with all sigma(i,i) = infinity (or a suitably large
       %number). 
       sigma_diag(i) = 1e7;
   end
end

sigma = sparse(1:numRxns, 1:numRxns, sigma_diag);

end

