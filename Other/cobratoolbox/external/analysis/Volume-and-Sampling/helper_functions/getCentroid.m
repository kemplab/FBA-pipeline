function [ centroid ] = getCentroid( P, eps)
%GETCENTROID compute an approximation of the centroid of P

%Input:

%P: a rounded polytope
%eps: error parameter

%Output:

%centroid: an estimate for the centroid

if nargin < 2
    eps = 0.1;
end

dim = size(P.A,2);
num_samples = round(dim/eps^2);
num_steps = dim+10;

points = genSamples(P,num_steps, num_samples);

centroid = mean(points,2);

end

