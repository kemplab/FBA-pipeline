function [ points ] = genSamplesExp( P,lambda, numSteps, numSamples, warmup )
%GENSAMPLES Compute samples from a polytope, taking numSteps steps of the
%Markov chain before saving a sample point

if nargin < 5
   warmup = 1e4; 
end

if nargin < 4
   numSamples = 1e3; 
end

if nargin < 3
    numSteps = 1e2;
end
%rotate so lambda is a coordinate vector
dim = size(P.A,2);
[q,~] = qr([lambda/norm(lambda)'; zeros(dim-1,1) eye(dim-1)]');
P.A = P.A*q;
P.N = P.N*q;
K = ConvexBody(P,[],.2,'');
x = zeros(K.dim,1);
resetSlacks(K,x);
points = zeros(size(P.N,1),numSamples);
lambda_rot = norm(lambda)*[1; zeros(dim-1,1)];
for i=1:warmup
   x = getNextPoint_Exp(K,x,lambda_rot,1); 
end

h = waitbar(0,'Computing samples...');
for i=1:numSamples
    for j=1:numSteps
        x = getNextPoint_Exp(K,x,lambda_rot,1);
    end
    points(:,i) = P.N*x+P.p_shift;
    waitbar(i/numSamples);
end
close(h);
end