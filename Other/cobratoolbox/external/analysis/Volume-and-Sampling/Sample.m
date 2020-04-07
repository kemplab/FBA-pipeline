function [x,T,steps] = Sample(P,E,eps,flags,start_pt)
%This function is a randomized algorithm to approximately sample from a convex
%body K = P \cap E with error parameter eps. The last 4 parameters are optional;
%you can see the default values below.

%---INPUT VALUES---
%P: the polytope, {x | P.A*x <= P.b}
%E: the ellipsoid, {x | (x-E.v)'E.E(x-E.v)<=1}
%eps: the target relative error
%flags: a string of input flags, see parseFlags.m
%start_pt: if you want a starting point for the walk

%---RETURN VALUES---
%x: the point approximately from the target distribution (default uniform)
%T: the rounding matrix. If no rounding, then T is identity matrix
%steps: the number of steps it took to reach the target point

%assign default values if not assigned in function call

if exist('flags','var')==0
    flags = '';
end
if exist('eps','var')==0
    eps = 0.20;
end
if exist('E','var')==0
    E=[];
end

%prepare our objects
K = ConvexBody(P,E,eps,flags);

%we don't need a starting Gaussian, but the following function computes the
%smallest enclosed ball
[~] = getStartingGaussian(K,0.1);

if exist('start_pt','var')==0
    %no point provided, pick one
    start_pt = zeros(K.dim,1);
elseif ~isempty(P)
    %we shifted the body to contain the origin, so also shift the start
    %point
    start_pt = start_pt - P.p;
end

if K.verb>=1
    fprintf('--------%d-Dimension Convex Body------\n\n', K.dim);
end

%let's initialize some helpful constants
[ratio,num_threads,C] = assignConstants(K);

if isKey(K.flagmap,'round')
    %rounding phase
    
    %round the body once as a preprocessing step
    %note that K is modified inside round()
    [T]=round(K,num_threads);
else
    %we are not rounding the body
    T=eye(K.dim);
end

if K.verb>=1
    fprintf('------Sample Start------\n');
end

%these are our starting points for each thread
%we shifted our body so that K contains the origin, so x \in K
x = start_pt;

if isKey(K.flagmap, 'a_stop')
    a_target = K.flagmap('a_stop');
else
    a_target = 0;
end

if isKey(K.flagmap, 'c_test')
    conv_test = K.flagmap('c_test');
else
    conv_test = 'sliding_window';
end

if isKey(K.flagmap,'plot_c')
   c_plotHandle = figure;
end

SCT = newSCT(conv_test, K.dim, eps, isKey(K.flagmap, 'plot_c'));

if isKey(K.flagmap, 'min_st')
    min_steps = K.flagmap('min_st');
else
    min_steps = 0;
end

if SCT.plotting==1
    c_pts = zeros(1000,1);
  else
        c_pts = [];
end

resetSlacks(K,x);

%the main loop. keep taking a step of our Markov chain until we're
%converged (by some statistical test)
while ~SCT.converged || SCT.steps < min_steps
    x = getNextPoint(K,x,a_target,1);
    
    SCT = updateSCT(SCT, conv_test, eps, x);
    
    if SCT.plotting==1
        if SCT.steps > size(c_pts,1)
           c_pts = [c_pts; zeros(size(c_pts,1),1)];
        end
        c_pts(SCT.steps) = SCT.plot_data;
    end
end

if K.verb>=2
    fprintf('Steps taken: %d', SCT.steps);
end

if K.verb>=1
    fprintf('------Sample End------\n');
end

steps = SCT.steps;

if SCT.plotting
    figure(c_plotHandle);
    plot(c_pts(1:SCT.steps));
    drawnow;
end
x=T*x;
end

function [ratio,num_threads,C] = assignConstants(K)
%initialize some hard-coded constants
if isKey(K.flagmap,'ratio')
    ratio = K.flagmap('ratio');
else
    ratio = 1-1/K.dim;
end
if isKey(K.flagmap,'num_t')
    num_threads = K.flagmap('num_t');
else
    num_threads = 5;
end
if isKey(K.flagmap,'C')
    C = K.flagmap('C');
else
    C = 2;
end

end

function [ SCT ] = newSCT(type, dim, curr_eps, to_plot)
%create the convergence test for generating a random sample

if strcmp(type,'sliding_window')==1
    %creating a new sliding window test
    SCT.W = ceil(4*dim^2+500);
    SCT.min_val=0;
    SCT.max_val=1;
    SCT.min_index=SCT.W;
    SCT.max_index=SCT.W;
    SCT.last_W = zeros(SCT.W,1);
    SCT.index = 1;
    %a random halfspace
    SCT.h = randn(dim,1);
    %the number which lie on one side of the halfspace
    SCT.h_num = 0;
    %number of steps we've taken
    SCT.steps = 0;
    %whether this test has converged
    SCT.converged = 0;
elseif strcmp(type,'mix')==1
    %creating a new test, mix for 8n^2 steps
    SCT.num_steps = 0;
    SCT.target_steps = 8*dim^2;
    SCT.converged = 0;
else
    error('Convergence test %s undefined.\n', type);
end

SCT.plotting = to_plot;

end

function [SCT] = updateSCT(SCT, type, eps, pt)

if strcmp(type,'sliding_window')==1
    SCT.h_num = SCT.h_num + (dot(pt,SCT.h)>0);
    SCT.steps = SCT.steps + 1;
    val = SCT.h_num/SCT.steps;
    SCT.last_W(SCT.index)=val;
    if val<=SCT.min_val
        SCT.min_val=val;
        SCT.min_index=SCT.index;
    elseif SCT.min_index==SCT.index
        [SCT.min_val,SCT.min_index]=min(SCT.last_W);
    end
    
    if val>=SCT.max_val
        SCT.max_val=val;
        SCT.max_index=SCT.index;
    elseif SCT.max_index==SCT.index
        [SCT.max_val,SCT.max_index]=max(SCT.last_W);
    end
    
    %check to see if the last W points are within sufficiently
    %small relative error. if so, stop this volume phase
%     SCT.plot_data = (SCT.max_val-SCT.min_val)/min(SCT.max_val,1-SCT.max_val);
    SCT.plot_data = val;
    if (SCT.max_val-SCT.min_val)/min(SCT.max_val,1-SCT.max_val)<=eps
        SCT.converged=1;
    end
    
    SCT.index = mod(SCT.index,SCT.W)+1;
elseif strcmp(type,'mix')==1
    SCT.num_steps = SCT.num_steps+1;
    if SCT.num_steps == SCT.target_steps
        SCT.converged = 1;
    end
    SCT.plot_data = SCT.num_steps;
else
    error('Convergence test %s undefined.\n', type);
end

end