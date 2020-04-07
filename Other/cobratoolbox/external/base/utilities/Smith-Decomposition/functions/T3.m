function T = T3(n,i,j,a)
T = eye(n);
if strcmp(class(a),'sym')
     T = sym(T);
end
    T(i,:) = T(i,:) + T(j,:)*a;
end