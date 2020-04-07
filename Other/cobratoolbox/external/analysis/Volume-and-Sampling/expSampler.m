function [samples,roundedP] = expSampler(P,numSkip,numSamples,options)
% EXPSAMPLER Generate uniform samples from the polytope P
%
%
%   polySampler will generate numSamples samples from P, taking
%   numSkip steps of a random walk between each sample
%
%   Rounding the polytope is a potentially expensive step. If you generate multiple rounds
%   of samples from a single model, you can save roundedP from the first round and
%   input it for subsequent rounds.
%
% INPUTS:
% P ... polytope object
% .A ...
% .b ... specifying Ax<=b
% .A_eq ...
% .b_eq ... specifying A_eq x <= b
% numSkip ... Number of steps of coordinate hit-and-run between samples
% numSamples ... Number of samples
%
% OPTIONAL INPUTS:
% options ...
% .toRound ... {0,1,2} Option to round the polytope before sampling.
%               0: no rounding
%               1 (default): round using max volume ellipsoid
%               2: round using isotropy (slower, more accurate)
% .roundedP ...{0,1} specifying this is a rounded polytope computed from a
%              previous run of the algorithm (default=0). If sampling from
%              the same model multiple times, it is highly recommended to
%              save roundedP and reuse it.
% .fullDim ... {0,1} option specifying if polytope is full-dimensional
%	       after restricting to P.A_eq = P.b_eq. Specifies whether
%	       we check for zero width facets in Ax<=b. default = 0
% .toPreprocess ... {0,1} option, default = 1, set to 0 if using
%                   roundedPolytope as input
% .bioModel ... {0,1} option stating if the input is a metabolic
%                     network. default = 0. allows for some optimizations
%                     if = 1 as these models have some special structure.
%
% OUTPUTS:
% samples ... n x numSamples matrix of random flux samples
% roundedPolytope ... The rounded polytope. Save for use in subsequent
%                     rounds of sampling.
%
% October 2016, Ben Cousins and Hulda S. Haraldsdóttir

% Define defaults
%if nargin==4 && exists
%    numSkip = 8*size(roundedP.A,2)^2;
%end

if nargin<3 || isempty(numSamples)
    numSamples = 1000;
end

if nargin<4
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
        [roundedP] = preprocessExp(P, options.toRound);
    else
        roundedP = P;
    end
end




%now we're ready to sample
samples = genSamplesExp(roundedP, numSkip, numSamples);

% samples = genSamplesGaussian(roundedPolytope,numSkip,numSamples,100*ones(size(roundedPolytope.N,1),1),eye(size(roundedPolytope.N,1)));

end




