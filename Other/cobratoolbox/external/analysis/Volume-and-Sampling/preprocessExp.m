function [P] = preprocessExp( P , lambda, options)
%PREPROCESS
% INPUTS:
% P ... polytope object
% .A ...
% .b ... specifying Ax<=b
% .A_eq ...
% .b_eq ... specifying A_eq x <= b
%
% OPTIONAL INPUTS:
% options ...
% .toRound ... {0,1,2} Option to round the polytope before sampling.
%               0: no rounding
%               1 (default): round using max volume ellipsoid
%               2: round using isotropy (slower, more accurate)
% .fullDim ... {0,1} option specifying if polytope is full-dimensional
%	       after restricting to P.A_eq = P.b_eq. Specifies whether
%	       we check for zero width facets in Ax<=b. default = 0
% .bioModel ... {0,1} option stating if the input is a metabolic
%                     network. default = 0. =1 allows for some minor
%                     optimizations (ignored if .fullDim=1)
%
% The algorithm will first check each inequality constraint to see if it has
% nonzero width. If it finds a width 0 facet, then the facet gets added
% to the equality subspace. The algorithm will then restrict to the
% equality subspace defined by P.A_eq, P.b_eq.
%
% After this preprocessing, the algorithms rounds the polytope for
% the volume/sampling algorithms to be accurate and efficient.
%
% Example usage:
%
% P = makeBody('long_box',10);
% P = preprocess(P);
% vol = Volume(P);

if nargin < 2
    error('Insufficient number of arguments');
end

if nargin<3
    options.t = 1;
end

if ~isfield(options,'toRound')
    options.toRound=1;
end
if ~isfield(options,'fullDim')
    options.fullDim=0;
end
if ~isfield(options,'bioModel')
    options.bioModel = 0;
end



dim = size(P.A,2);

%check for width 0 facets to make sure we are full dimensional
if options.fullDim==0
    
    fprintf('Checking for width 0 facets...\n');
    
    eps_cutoff = 1e-7;
    %check if we are doing a metabolic network
    if options.bioModel==1
        %this utilizes that the polytope is a subspace intersected with
        %a box
        if exist('fastFVA')==2
            %check if we can do parallel for fastFVA
            v=ver;
            PCT='Parallel Computing Toolbox';
            if  any(strcmp(PCT,{v.Name}))
                p = parcluster('local');
                setWorkerCount(p.NumWorkers);
            end
            [minFlux, maxFlux] = fastFVA(options.model,100);
        else
            [minFlux, maxFlux] = fluxVariability(options.model);
        end

        
        tol = 1e-6;
        optPercentage = 100;
        optSol = optimizeCbModel(options.model,'max', 0, true);
        objValue = floor(optSol.f/tol)*tol*optPercentage/100;
        P.A = [P.A; -options.model.c'];
        P.b = [P.b; -objValue];
        
        
        isEq = (maxFlux - minFlux) < eps_cutoff;
        eq_constraints = sparse(sum(isEq),size(P.A_eq,2));
        eq_constraints(:,isEq) = speye(sum(isEq));
        
        fprintf('Found %d degenerate reactions, adding them to the equality subspace.\n', sum(isEq));
        
        P.A_eq = [P.A_eq; eq_constraints];
        P.b_eq = [P.b_eq; minFlux(isEq)];
    else
        %check the width of every facet
        [widths, vals] = getWidths(P);
        
        
        num_eq = sum(abs(widths)<eps_cutoff);
        
        fprintf('Found %d width 0 facets, adding them to the equality subspace.\n', num_eq);
        if num_eq > 0
            
            eq_constraints = zeros(num_eq,dim);
            eq_rhs = zeros(num_eq,1);
            curr = 1;
            for i=1:length(widths)
                if abs(widths(i)) < eps_cutoff
                    %this facet has width 0
                    %we add this constraint to the equality
                    %subspace, with vals(i) defining
                    %the shift of the facet
                    
                    eq_constraints(curr,:) = P.A(i,:);
                    eq_rhs(curr) = vals(i);
                    curr = curr + 1;
                end
            end
            
            P.A_eq = [P.A_eq; eq_constraints];
            P.b_eq = [P.b_eq; eq_rhs];
        end
    end
    
end


%restrict to the degenerate subspace, if it exists
if isfield(P,'A_eq') && ~isempty(P.A_eq)
    %check if we can use better null space code
    if exist('getNullSpace')==2
        N = getNullSpace(P.A_eq,0);
    else
        %otherwise just use matlab's built-in one
        warning('Using MATLAB''s null(), results may be inaccurate. Recommend installing LuSOL.\n');
        N = null(P.A_eq);
    end
    
    %in case we want the volume dilation factor when restricting to the null space, uncomment
    %the below line (e.g. if the product of all the singular values of N = 1, then there is no
    %volume change)
    
    P.vol_increase = 1/prod(svd(full(N)));
    
    %find a point in the null space to define the shift
    z = P.A_eq \ P.b_eq;
    N_total = N;
    p_shift = z;
    P.b = P.b - P.A * z;
    P.A = P.A*N;
    dim = size(P.A,2);

    fprintf('Now in %d dimensions after restricting\n', dim);
else
    N_total = eye(dim);
    p_shift = zeros(dim,1);
    P.vol_increase = 1;
end

%remove zero rows
row_norms = sqrt(sum(P.A.^2,2));
P.A = P.A(row_norms>1e-6,:);
P.b = P.b(row_norms>1e-6,:);

%maybe add LP presolve here

fprintf('Removed %d zero rows\n', sum(row_norms>1e-6));

%scale the matrix to help with numerical precision
dim = size(P.A,2);
if exist('gmscale')==2
    fprintf('Preconditioning A with gmscale\n');
    [cs,rs] = gmscale(P.A,0,0.99);
    P.A = diag(1./rs)*P.A*diag(1./cs);
    P.b = diag(1./rs)*P.b;
    N_total = N_total*diag(1./cs);
end

row_norms = sqrt(sum(P.A.^2,2));
P.A = diag(1./row_norms)*P.A;
P.b = diag(1./row_norms)*P.b;

fprintf('Rounding...\n');

Q = P;
% tol = 1e-6;
% optPercentage = 50;
LP.A = P.A;
LP.b = P.b;
LP.c = lambda' * N_total;
LP.lb = -Inf * ones(length(LP.c),1);
LP.ub = Inf*ones(length(LP.c),1);
LP.osense = -1;
LP.csense = repmat('L',size(LP.A,1),1);
solution = solveCobraLP(LP);
Q.A = [Q.A; -LP.c];
Q.b = [Q.b; -solution.obj/2/norm(LP.c)];

if options.toRound==1
    
    T = eye(dim);
    if dim<200
        %we can just solve the interior point method once
        if exist('solveCobraLP')==2
            [x0,dist] = getCCcenter(Q.A,Q.b);
        else
            %let mve_run use matlab's lp solver to select a starting point
            [~,x0] = mve_presolve_cobra(Q.A,Q.b,150,1e-6);
        end
        
        [T_shift, Tmve,converged] = mve_run_cobra(Q.A,Q.b, x0,1e-8);
        
        [P,~, ~, ~] = shiftPolytope(P, N_total, p_shift, T, Tmve, T_shift);
        [Q,N_total, p_shift, T] = shiftPolytope(Q, N_total, p_shift, T, Tmve, T_shift);
        if converged~=1
            fprintf('There was a problem with finding the maximum volume ellipsoid.\n');
        end
    else
        %if dimension is large, be a little more careful
        %the below loop looks silly for an interior point method, but is actually
        %quite important for numerical stability. while normally you'd only call an optimization routine
        %once, we call it iteratively--where in each iteration, we map a large
        %ellipsoid to the unit ball. the idea is that the "iterative roundings"
        %will make each subsequent problem easier to solve numerically.
        max_its = 20;
        its = 0;
        reg=1e-3;
        Tmve = eye(dim);
        converged = 0;
        while (max(eig(Tmve))>6*min(eig(Tmve)) && converged~=1) || reg>1e-6 || converged==2
            tic;
            its = its+1;
            %check if we can use the good lp solver
            if exist('solveCobraLP')==2
                [x0,dist] = getCCcenter(Q.A,Q.b);
            else
                %let mve_run use matlab's lp solver to select a starting point
                [~,x0] = mve_presolve_cobra(Q.A,Q.b,150,1e-6);
            end
            
            reg = max(reg/10,1e-10);
            [T_shift, Tmve,converged] = mve_run_cobra(Q.A,Q.b, x0,reg);
            [P,~, ~, T] = shiftPolytope(P, N_total, p_shift, T, Tmve, T_shift);
            [Q,N_total, p_shift, T] = shiftPolytope(Q, N_total, p_shift, T, Tmve, T_shift);
            row_norms = sqrt(sum(Q.A.^2,2));
            Q.A = diag(1./row_norms)*Q.A;
            Q.b = diag(1./row_norms)*Q.b;
            if its==max_its
                break;
            end
            
            fprintf('Iteration %d: reg=%.1e, ellipsoid vol=%.1e, longest axis=%.1e, shortest axis=%.1e, x0 dist to bdry=%.1e, time=%.1e seconds\n', its, reg, det(Tmve), max(eig(Tmve)), min(eig(Tmve)), dist, toc);
        end
        
        if its==max_its
            fprintf('Reached the maximum number of iterations, rounding may not be ideal.\n');
        end
    end
    
    
    if min(P.b)<=0
        if exist('solveCobraLP')==2
            [x,~] = getCCcenter(P.A,P.b);
        else
            %let mve_run use matlab's lp solver to select a starting point
            [~,x0] = mve_presolve_cobra(P.A,P.b,150,1e-6);
        end
        [P,N_total,p_shift,T] = shiftPolytope(P,N_total,p_shift,T,eye(dim),x);
        fprintf('Shifting so the origin is inside the polytope...rounding may not be ideal.\n');
    else
        fprintf('Maximum volume ellipsoid found, and the origin is inside the transformed polytope.\n');
    end
    
    P.p = zeros(dim,1);
    
elseif options.toRound==2
    K = ConvexBody(P,[],.1,'');
    T = round(K,5);
    P.p = zeros(dim,1);
else
    
    if exist('solveCobraLP')==2
        [P.p,~] = getCCcenter(P.A,P.b);
    else
        %let mve_run use matlab's lp solver to select a starting point
        [~,P.p] = mve_presolve_cobra(P.A,P.b,150,1e-6);
    end
    T = eye(size(N_total,2));
end

P.N = N_total;
P.p_shift = p_shift;
P.T = T;
P.vol_increase = P.vol_increase / abs(det(P.T));
P.lambda = -lambda' * P.N;
end

%compute the center of the Chebyshev ball in the polytope Ax<=b
function [CC_center,radius] = getCCcenter(A,b)
dim = size(A,2);
a_norms = sqrt(sum(A.^2,2));

LP.A = [A a_norms];
LP.b = b;
LP.c = [zeros(dim,1); 1];
LP.lb = -Inf * ones(dim+1,1);
LP.ub = Inf*ones(dim+1,1);
LP.osense = -1;
LP.csense = repmat('L',size(LP.A,1),1);
solution = solveCobraLP(LP);
if solution.stat == 1
    CC_center = solution.full(1:dim);
    radius = solution.obj;
else
    solution
    error('Could not solve the LP, consult the information above.');
end
end

%shift the polytope by a point and apply a transformation, while retaining
%the information to undo the transformation later (to recover the samples)

%let x denote the original space, y the current space, and z the new space
%we have
%
%   P.A y <= P.b   and x = N*y+p
%
%  applying the transformation
%
%   trans * z + shift = y
%
% yields the system
%
%  x = (N * trans) * z + N*shift + p
%  (P.A * trans) * z <= P.b - P.A * shift
function [P,N,p,T] = shiftPolytope(P,N,p,T,trans, shift)
p = p + N*shift;
N = N * trans;
T = T * trans;
P.b = P.b - P.A*shift;
P.A = P.A*trans;
end

function [widths, vals] = getWidths(P)

options = optimset('Display','none');
warning('off');

widths = zeros(size(P.A,1),1);
vals = zeros(size(P.A,1),1);

fprintf('%d constraints to check.\n', length(widths));
chunk = ceil(length(widths)/10);

for i=1:length(widths)
    %disp(i)
    [x,max_dist] = linprog(-P.A(i,:), P.A,P.b,P.A_eq,P.b_eq,[],[],[],options);
    [~,min_dist] = linprog(P.A(i,:), P.A,P.b,P.A_eq,P.b_eq,[],[],[],options);
    
    if mod(i,chunk)==0
        fprintf('%d/%d done..\n', i, length(widths));
    end
    
    widths(i) = abs(max_dist+min_dist)/norm(P.A(i,:),2);
    vals(i) = P.A(i,:)*x;
    
end

end