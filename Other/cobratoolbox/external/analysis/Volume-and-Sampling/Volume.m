function [volume,T,steps] = Volume(P,E,eps,flags)
%This function is a randomized algorithm to approximate the volume of a convex
%body K = P \cap E with relative error eps. The last 3 parameters are optional; 
%if not provided, default values are assigned (see below).

%---INPUT VALUES---
%P: the polytope, {x | P.A*x <= P.b}
%E: the ellipsoid, {x | (x-E.v)'E.E(x-E.v)<=1}
%eps: the target relative error
%flags: a string of input flags, see parseFlags.m

%---RETURN VALUES---
%volume: the computed volume estimate
%T: the rounding matrix. If no rounding, then T is identity matrix
%steps: the number of steps the volume algorithm took

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
if K.verb>=1
    fprintf('--------%d-Dimension Convex Body------\n\n', K.dim);
end

if isKey(K.flagmap,'plot')
    plotting = 1;
    plotHandle = figure;
else
    plotting = 0;
end

if isKey(K.flagmap,'plot_c')
   c_plotHandle = figure;
end

%let's initialize some helpful constants
[ratio,num_threads,C] = assignConstants(K);

if isKey(K.flagmap,'round')
    %rounding phase  
    
    if strcmp(K.flagmap('round'),'mve')==1
        
        fprintf('here\n');
        if ~isempty(E)
            error('Cannot round with MVE if there is an ellipsoidal constraint.');
        end
        options.fullDim = 1;
        options.toRound=1;
        Q.A = P.A;
        Q.b = P.b;
        roundedP = preprocess(Q,options);
        K.A = roundedP.A;
        K.b = roundedP.b;
        K.p = zeros(K.dim,1);
        T = roundedP.T;
    else
        %note that K is modified inside round
        [T] = round(K,num_threads);
    end
else
    %we are not rounding the body
    T=eye(K.dim);
end

if K.verb>=1
    fprintf('------Volume Start------\n');
end

%compute the annealing schedule that keeps E(Y^2)/E(Y)^2<=C
[a_sched] = getAnnealingSchedule(K,ratio,num_threads,C);
K.m = length(a_sched);

%these are our starting points for each thread
%we shifted our body so that K contains the origin, so x \in K
x = zeros(K.dim,num_threads);

resetSlacks(K,x);

%compute the initial volume, multiplied by the determinant of the rounding
% matrix.
volume = (pi/a_sched(1))^(K.dim/2)*abs(det(T));

%check if we should increase the volume, in case the preprocessing
%algorithm distorted it
if isfield(P,'vol_increase')
    volume = volume*P.vol_increase;
end

%initialize additional helpful variables
fn = zeros(length(a_sched),1);
its = zeros(length(a_sched),1);

if K.verb>=1
    fprintf('Num Phases: %d\n', K.m);
end

if isKey(K.flagmap, 'c_test')
    conv_test = K.flagmap('c_test');
else
    %default to sliding window convergence test
    conv_test = 'sliding_window';
end

for i=1:length(a_sched)-1
    if K.verb>=1
        fprintf('Phase %d Volume: %e', i-1,volume);
        if K.verb>=1
            %compute how many "ratios" we've stepped down so far
            sigma_index=round(log(a_sched(i)/a_sched(1))/log(ratio));
            fprintf(',   sigma_%d=%e', sigma_index,1/sqrt(2*a_sched(i)));
        end
        fprintf('\n');
    end
    %prepare to sample this volume phases
    
    %allocate the error to this phase
    curr_eps = K.eps/sqrt(K.m);
    
    %initialize our convergence test
    VCT = newVCT(conv_test, K.dim, curr_eps, isKey(K.flagmap, 'plot_c'));
    
    %check if we should store data for plotting
    if plotting==1
       %we just store the first 2 dimensions, since that's all we will
       %be displaying
       pts = zeros(1000,2);
    else
        pts = [];
    end
    
    if VCT.plotting==1
        c_pts = zeros(1000,1);
    else
        c_pts = [];
    end
    
    %this is the flag for when we converge
    done=0;
    if isKey(K.flagmap, 'min_st')
       min_steps = K.flagmap('min_st');
    else
        min_steps = 0;
    end
    
    %keep taking steps until we converge
    while ~done || its(i) < min_steps
        %advance each of the threads one step
        for j=1:num_threads
            
            %take one step of our Markov chain
            x(:,j) = getNextPoint(K,x(:,j),a_sched(i),j);
            
            if ~in_K(K,x(:,j))
               error('Found a point not in K.'); 
            end
            
            %update the current volume ratio with this new point
            its(i) = its(i)+1;
            fn(i) = fn(i) + eval_exp(x(:,j),a_sched(i+1))/eval_exp(x(:,j),a_sched(i));
            
            %Update the sliding window, keeping track of min/max over last W
            %steps. For random points, this will work in O(1) expected amortized time.
            val = fn(i)/its(i);
            [VCT] = updateVCT(VCT, conv_test, curr_eps, val, x(:,j));
            
            %store this point if we're plotting the results
            if plotting==1
                if its(i)>size(pts,1)
                   pts = [pts; zeros(size(pts,1), 2)];
%                    num_pts = 2*num_pts; 
                end
                pts(its(i),1:2) = x(1:2,j)';
            end
            
            %check if we're plotting the convergence test
            
            if VCT.plotting==1
               if its(i)>size(c_pts,1)
                   c_pts = [c_pts; zeros(size(c_pts,1),1)];
               end
               c_pts(its(i)) = VCT.plot_data;
            end
            
            if VCT.converged
               done = 1;
               break; 
            end
            
        end
        
    end
    if plotting == 1
        PlotDistribution(pts(1:its(i),:),plotHandle);
    end
    
    if VCT.plotting == 1
       figure(c_plotHandle);
       plot(c_pts(1:its(i))); 
       drawnow;
    end
    if K.verb>=2
        fprintf('Steps in Phase %d: %d\n', i-1, its(i));
    end
    volume = volume*fn(i)/its(i);
end

if K.verb>=1
    fprintf('------Volume Complete------\n\n');
    fprintf('Final Volume: %e,  final sigma=%e\n',volume, 1/sqrt(2*a_sched(end)));
    fprintf('Total steps: %d\n', sum(its));
end

steps=sum(its);
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

function [ VCT ] = newVCT(type, dim, curr_eps, to_plot)
%create a convergence test, which monitors the volume and decides
%if it has converged in this phase

if strcmp(type,'sliding_window')==1
      %creating a new sliding window test
       VCT.W = ceil(4*dim^2+500);
        VCT.min_val=-1e100;
        VCT.max_val=1e100;
        VCT.min_index=VCT.W;
        VCT.max_index=VCT.W;
        VCT.last_W = zeros(VCT.W,1);
        VCT.index = 1;
        VCT.converged = 0;
else
    error('Convergence test %s undefined.\n', type);
end

VCT.plotting = to_plot;

end

function [VCT] = updateVCT(VCT, type, curr_eps, val, pt)
%Update the volume convergence test after a step of the random walk.

    if strcmp(type,'sliding_window')==1
        
         VCT.last_W(VCT.index)=val;
            if val<=VCT.min_val
                VCT.min_val=val;
                VCT.min_index=VCT.index;
            elseif VCT.min_index==VCT.index
                [VCT.min_val,VCT.min_index]=min(VCT.last_W);
            end
            
            if val>=VCT.max_val
                VCT.max_val=val;
                VCT.max_index=VCT.index;
            elseif VCT.max_index==VCT.index
                [VCT.max_val,VCT.max_index]=max(VCT.last_W);
            end
            
            %check to see if the last W points are within sufficiently
            %small relative error. if so, stop this volume phase
%             VCT.plot_data = (VCT.max_val-VCT.min_val)/VCT.max_val;
            VCT.plot_data = val;
            if (VCT.max_val-VCT.min_val)/VCT.max_val<=curr_eps/2
                VCT.converged=1;
            end
            
            VCT.index = mod(VCT.index,VCT.W)+1;
    else
       error('Convergence test %s undefined.\n', type); 
    end
            
end