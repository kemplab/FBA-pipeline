function [ samples ] = gaussianSampler( P, numSkip, numSamples, mu, sigma, options )
%GAUSSIANSAMPLER Summary of this function goes here
%   Detailed explanation goes here

if nargin<3 || isempty(numSamples)
    numSamples = 1000;
end

if nargin<6
   options = []; 
end

if ~isfield(options,'toRound')
    options.toRound=1;
end

if ~isfield(options,'roundedP')
    options.roundedP=0;
end
if ~isfield(options,'fullDim')
    options.fullDim=0;
end
if ~isfield(options,'toPreprocess')
    options.toPreprocess=1;
end
if ~isfield(options,'bioModel')
    options.bioModel = 0;
end


% Preprocess model
if options.toPreprocess==1
    
    %check to make sure P.A and P.b are defined, and appropriately sized
    if (isfield(P,'A')==0 || isfield(P,'b')==0) || (isempty(P.A) || isempty(P.b))
        %either P.A or P.b do not exist
        error('You need to define both P.A and P.b for a polytope {x | P.A*x <= P.b}.');
    end
    
    [num_constraints,dim] = size(P.A);
    
    if exist('numSkip')~=1 || isempty(numSkip)
        numSkip=8*dim^2;
    end
    
    fprintf('Currently (P.A, P.b) are in %d dimensions\n', dim);
    
    if size(P.b,2)~= 1 || num_constraints ~= size(P.b,1)
        error('Dimensions of P.b do not align with P.A.\nP.b should be a %d x 1 vector.',num_constraints);
    end
    if (isfield(P,'A_eq')==0 || isempty(P.A_eq)) && ...
            (isfield(P,'b_eq')==0 || isempty(P.b_eq))
        P.A_eq = [];
        P.b_eq = [];
    end
    
    %preprocess the polytope of feasible solutions
    %restict to null space
    %round via maximum volume ellipsoid
    %     save('recon_widthed.mat', 'P');
    if options.roundedP==0
        [roundedP] = preprocess(P, options.toRound);
    else
        roundedP = P;
    end
end


%now we're ready to sample
% samples = genSamples(roundedP, numSkip, numSamples);

samples = genSamplesGaussian(roundedP,numSkip,numSamples,mu, sigma);

end

