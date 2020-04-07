function [ dist ] = ellipsoidDist(Q, y)
global mat pt
mat = Q;
pt = y;
%ELLIPSOIDDIST compute the minimum distance from the point y
% to a point x such that x^T Q x = 1.
% We use the method of Lagrange multipliers.
% The optimization problem is
%    min d(x,y)^2
%    s.t.
%      x^T Q x = 1
%
%The Langrangian is L(x,lambda) = ||x-y||^2 + lambda * (x^T Q x - 1).
%The partial derivatives of L are
%  L_x(x,lambda) = 2*(x-y) + lambda*(Q+Q^T)x
%  L_lambda(x,lambda) = x^T Q x - 1
%
%We need all these partial derivatives equal to zero. We use MATLAB's 
%multivariate polynomial root finder to find a solution (hopefully).
%

%First, build our linear system.
%
% lower = 0;
% upper = 1e50;
%
% while upper-lower>1e-3
%    mid = (upper+lower)/2;
%
%    if eval_fun(mid,Q,y)<=1
%       %guess was too high
%       fprintf('Guess %e too high, eval_fun = %e\n', mid, eval_fun(mid,Q,y));
%       upper = mid;
%    else
%       fprintf('Guess %e too low, eval_fun = %e\n', mid, eval_fun(mid,Q,y));
%        lower = mid;
%    end
% end
%
% % dim = length(y);
% %
% % A = (Q + Q' + 2*eye(dim));
% %
% % b = 2*y;
% %
% % x = linsolve(A,b);
% %
% % lambda = sqrt(x' * Q * x);
% %
% % x = x / lambda;
% %
% % dist = norm(x-y);
%
%
% end
%
%
% %evaluate
% function [res] = eval_fun(lambda,Q,y)
%     A = lambda*(Q + Q') + 2*eye(length(y));
%     b = 2*y;
%
%     x = linsolve(A,b);
%
%     res =  x'*Q*x;
%
% end

yT_Q_y = y' * Q * y;

scale = sqrt(yT_Q_y);

if scale ~= 0
    x0 = y / scale;
    
    x0 = [x0; 0];
    options = optimoptions('fsolve','Display','off');
    [x,fval] = fsolve(@myfun,x0,options);
    
%     x = x(1:end-1);
    dist = norm(x(1:end-1)-y);
else
    eigen_vals=eig(Q);
    dist = 1/sqrt(max(eigen_vals));
end

end

function F = myfun(x)
global mat pt
lambda = x(end);
z = x(1:end-1);

F = zeros(length(x),1);
F(1:end-1) = lambda * (mat+mat') * z + 2*(z-pt);
F(end) = z' * mat * z -1;
end