function [T,D,S] = smithFormInt(A)
% Smith normal form of an integer matrix.
iters = max(size(A))-1;
T = eye(size(A,1));
S = eye(size(A,2));
D = A;
i = 1;

% Apply row and column operations 
while i<iters+1
    [T LD] = findgcdReduceCol(D,i-1,T);
    [D S] = findgcdReduceRow(LD,i-1,S);
    i = i + 1;
    %if row operation messed up column repeat step
    Col_check =  D(i:end,i-1);
    if any(Col_check(:)) && ~isempty(Col_check)
        i = i-1;
    end
end
%Check last diagonal element
if D(end,end) < 0
    T = T1(size(D,1),size(D,1),-1)*T;
    D = T1(size(D,1),size(D,1),-1)*D;
end

% Decomposition Check
Check = T*A*S - D;
if any(Check(:))
    error('Something is wrong!!!!!!!!@%@$%$#^')
end
end