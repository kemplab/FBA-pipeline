function T = T2(n,i,j)
T = eye(n);
T([i j],:) = T([j i],:);
end