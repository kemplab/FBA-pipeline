function T = T1(n,i,c)
T = eye(n);
if strcmp(class(c),'sym')
     T = sym(T);
end
T(i,:) = T(i,:)*c;
end