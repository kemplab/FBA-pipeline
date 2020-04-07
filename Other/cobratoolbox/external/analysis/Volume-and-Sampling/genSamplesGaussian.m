function [ points, P ] = genSamplesGaussian( P,mu, sigma2, numSteps, numSamples, warmup)
%GENSAMPLES Compute samples from a polytope, taking numSteps steps of the
%Markov chain before saving a sample point

if nargin < 6
    warmup = 1e4;
end

if nargin < 5
   numSamples = 1e3; 
end

if nargin < 4
    numSteps = 1e2;
end


%Let y denote projected space, x original.
%We have x = P.N * y + P.p_shift, and the polytope corresponding to y has 
%maximum volume ellipsoid = unit ball.
%We want to restrict the Gaussian to the projected space
%To do this, we need B,z such that y = B * x + z
%Since P.N is not square, we use the SVD to get the "inverse" of P.N

[U,S,V] = svd(P.N);
s = diag(S);
Sinv = diag(1./s);
Sinv = [Sinv zeros(size(V,1),size(U,1)-size(V,1))];

%compute B s.t. B*P.N = I
%so now y = B*x + z
B = V*Sinv*U';
z = -B * P.p_shift;

%generate the new Gaussian parameters using mu, sigma2 and the
%transformation matrix and shift
new_mu = B * mu + z;
new_sigma2 = P.N' * sigma2 * P.N;
%look at the variances
[V,D] = eig(new_sigma2); %new_sigma2 = V * D * V'

%Apply a transformation so that the variance in each direction is at least
%1/2 (where 1/2 is a somewhat arbitrary constant below 1).
%The large variance directions will already be bounded due to the fact
%that the Max Vol Ellipsoid in the set is the unit ball.
%We apply this transformation since the variance could be arbitrarily small
%in a direction, which can subsequently make the sampling take arbitrarily
%long to converge.

%Let x denote the original variable, and let y denote the transformed space
%Going to transform y = V * F * V' * x where V is computed via the eigen
%decomposition of new_sigma2, and F is a diagonal matrix that rescales the
%variances
% F = eye(size(D,1));

F = zeros(size(D));
var_threshold = 1/2;
for i=1:size(D,1)
    if D(i,i) < var_threshold
%         F(i,i) = 1/2;
        F(i,i) = sqrt(D(i,i)/var_threshold);
    else
        F(i,i) = D(i,i);
    end
end

% [V,D] = eig(rounded_sigma);


%shift so that the Gaussian is centered at 0 (makes future calculations
%easier) y-mu -> y1
P.b = P.b - P.A*new_mu;
P.p_shift = P.p_shift + P.N * new_mu;
P.p = -new_mu;
new_mu = 0*new_mu;

%now transform using y = V*F*V'*y'
P.A = P.A * V * F * V';
P.N = P.N * V * F * V';
new_mu = V*diag(1./diag(F))*V'*new_mu;
% P.p_shift = V*diag(1./diag(F))*V'*P.p_shift;

%Let's rotate the space again so that the covariance matrix is diagonal.
%This is useful to the later sampling procedure so that when we ask for the
%variance along an axis-aligned chord, it is a constant time operation.

%(maybe come back later and merge the below step with the above, but
%keeping the code verbose for now for readability)
rounded_sigma = diag(1./diag(F))*D*diag(1./diag(F));
% rounded_sigma = F;
P.A = P.A*V;
P.N = P.N*V;
new_mu = V'*new_mu;
% P.p_shift= V'*P.p_shift;

%a strictly interior point
%this is the point in the center of the ball contained in P
P.p = -new_mu;

%find the point in our convex body that has the highest Gaussian density
%this should be a reasonable starting point for a random walk
K = ConvexBody(P,[],.2,'');
QPproblem.A = P.A;
QPproblem.b = P.b;
QPproblem.csense = 'L';
QPproblem.osense = 1;
dim = size(P.A,2);
QPproblem.lb = -Inf*ones(dim,1);
QPproblem.ub = Inf*ones(dim,1);
QPproblem.c = zeros(dim,1);
QPproblem.F = diag(1./diag(rounded_sigma));
solution = solveCobraQP(QPproblem);

% x = solution.full;

x = (dim-1)/dim * solution.full + 1/dim * P.p; %give a little padding to make the current point strictly inside the polytope
% x = P.p;
resetSlacks(K,x);
points = zeros(size(P.N,1),numSamples);

sigma_v = diag(rounded_sigma);
min(P.A*x<P.b)
for i=1:warmup
   x = getNextPoint_MVnormal(K,x,sigma_v,1); 
end

for i=1:numSamples
    if mod(i,round(numSamples/20))==0
        fprintf('%d/%d\n', i, numSamples);
    end
    for j=1:numSteps
        [x] = getNextPoint_MVnormal(K,x,sigma_v,1);
    end
    
    %make sure P.A*x <= P.b
%     if min(P.A*x<=P.b)==0
%         fprintf('uh-oh\n');
%     end
    
    points(:,i) = P.N*x+P.p_shift;
end

end