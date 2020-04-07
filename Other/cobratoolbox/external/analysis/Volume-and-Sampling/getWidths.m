function [widths, vals, minFlux, maxFlux] = getWidths(P)

options = optimset('Display','none');

widths = zeros(size(P.A,1),1);
vals = zeros(size(P.A,1),1);

minFlux = zeros(size(widths));
maxFlux = zeros(size(widths));

fprintf('%d constraints to check.\n', length(widths));
chunk = ceil(length(widths)/10);

% LP.A = [P.A; P.A_eq];
% LP.b = [P.b+1e-8; P.b_eq];
% LP.csense = [repmat('L',size(P.b)); repmat('E',size(P.b_eq))];
% LP.osense = 1;
% LP.lb = -Inf*ones(size(P.A,2),1);
% LP.ub = -LP.lb;
LP.A = sparse(P.A_eq);
LP.b = P.b_eq;
LP.csense = repmat('E',size(LP.A,1),1);
LP.osense = 1;
LP.ub = P.b(1:size(P.A,2));
LP.lb = -P.b(size(P.A,2)+1:end);
%LP.basis = [];
LP.c = zeros(size(LP.A,2),1);
LP.c(1)=1;
parameters.FeasibilityTol = 1e-8;

% tempSolution = gurobi(LP);
%tempSolution = solveCobraLP(LP);
%LP.basis = tempSolution.basis;
% LP.vbasis = tempSolution.vbasis;
% LP.cbasis = tempSolution.cbasis;

% LP_pre = gurobi_read('presolved_recon.mps');

for i=1:length(widths)
    %disp(i)
    %[~,max_dist] = linprog(f, [A1; K.A(i,:) -1],[K.b; K.b(i)],[],[],[],[],[],options);
    %[~,min_dist] = linprog(f, [A1; K.A(i,:) 1],[K.b; K.b(i)],[],[],[],[],[],options);
    
    if mod(i,chunk)==0
        fprintf('%d/%d done..\n', i, length(widths));
    end
    
    %     tic;
    %     [x,max_dist,exitflag] = linprog(-P.A(i,:), P.A, P.b,[],[],[],[],[],options);
    %     toc
    LP.c = full(P.A(i,:));
%     solution = gurobi(LP_pre);
    [solution] = solveCobraLP(LP, parameters);
    
    %LP.basis = solution.basis;
    
    x = solution.full;
    max_dist = solution.obj;
    
    %     tic;
    %     [x,~,~,~,max_dist,~,~,~,basis] = gurobi(-P.A(i,:), P.A, lb, ub, b_L, P.b);
    %     toc
    
    if solution.stat == 0
        error('The polytope provided has no feasible points! (so it appears)');
    end
    %     [y,min_dist] = linprog(P.A(i,:), P.A, P.b,[],[],[],[],[],options);
    
    LP.c = -LP.c;    
%     solution = gurobi(LP_pre);
    [solution] = solveCobraLP(LP, parameters); %,'feasTol',1e-6,'optTol',1e-6);
    

    %LP.basis = solution.basis;
    
    min_dist = solution.obj;
    
    if solution.stat == 0 || isempty(min_dist)
        error('The polytope provided has no feasible points! (so it appears)');
    end
    
    widths(i) = abs(max_dist+min_dist)/norm(P.A(i,:),2);
    vals(i) = P.A(i,:)*x;
    %     if widths(i) < 1e-7
    %         fprintf('(%e,%e)\n', vals(i), min_dist);
    %     end
end

end