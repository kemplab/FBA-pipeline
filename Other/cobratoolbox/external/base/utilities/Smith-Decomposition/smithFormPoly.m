function [T,D,S] = smithFormPoly(A)
% Smith normal form of a matrix with polynomial entries
iters = max(size(A))-1;
T = eye(size(A,1));
S = eye(size(A,2));
D = A;
i = 1;

% Apply row and column operations 
while i<iters+1
    [T LD] = findgcdReduceColPoly(D,i-1,T);
    [D S] = findgcdReduceRowPoly(LD,i-1,S);
    i = i + 1;
    %if row operation messed up column repeat step
    Col_check =  sym2poly(D(i:end,i-1));
    if any(Col_check(:)) && ~isempty(Col_check)
        i = i-1;
    end
%     break;
end

% Check if last poly is monic
inv_factors = diag(D);
last = sym2poly(inv_factors(end));
if last(1)~=1 && size(A,1)>1 && last(1)~=0
    T = T1(size(A,1),length(inv_factors),1/last(1))*T;
    D = T1(size(A,1),length(inv_factors),1/last(1))*D;
end

TAS = T*A*S;
TAS = precisionCheck(TAS);

% Decomposition Check
if ((TAS - D) ~= 0)
    error('Something is wrong!!!!!!!!@%@$%$#^')
end

end