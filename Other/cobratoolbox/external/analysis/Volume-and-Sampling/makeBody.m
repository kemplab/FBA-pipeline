function [body, volume] = makeBody( shape, dim )
%MAKEBODY This function will create some example bodies that can be used to
%call "Volume". "body" will be the body, either a polytope or ellipsoid.

%A polytope is given in the form [A b], {x | Ax <= b}. The polytope will be
%an m x (n+1) matrix, for m constraints and n dimensions.
%An ellipsoid is given in the form [E v], {x | (x-v)'E(x-v)<=1}. The
%ellipsoid will be an n x (n+1) matrix, for n dimensions.

if strcmp(shape,'ball')==1
    body.E = eye(dim);
    body.v = zeros(dim,1);
    volume=pi^(dim/2)/gamma(dim/2+1);
elseif strcmp(shape,'ellipsoid')==1
    body.E = eye(dim);
    body.E(dim,dim)=1e-4;
    body.v = zeros(dim,1);
    volume=pi^(dim/2)/gamma(dim/2+1)*100;
elseif strcmp(shape,'birkhoff')==1
    %describe system as Ax<=b and Cx=d
    %need to get subspace spanned by Cd 
    %and apply that transformation to Ax<=b
    %b starts out as all zeros, d as all ones
    A=-eye(dim^2);
    C=zeros(2*dim,dim^2);
    for i=1:dim
        %row constraints
        for j=(i-1)*dim+1:i*dim
          C(i,j)=1; 
        end
        
        %col constraints
        for j=i:dim:dim^2
           C(dim+i,j)=1; 
        end
    end

    nullC=null(C);
    %let v be any solution to Cx=d
    v=ones(dim^2,1)/dim;
    
    new_A=A*nullC;
    new_b=zeros(dim^2,1)-A*v;
body.A = new_A;
body.b = new_b;
    body.p=zeros((dim-1)^2,1);
    volume = -1;
elseif strcmp(shape,'triply_stochastic')==1
    A = -eye(dim^3);
    C = zeros(3*dim, dim^3);
    for i=1:dim
       for j=(i-1)*dim^2+1:i*dim^2
           C(i,j)=1;
       end
       
       for j=i:dim:dim^3
          C(dim+i,j)=1; 
       end
       
       for j=1:dim
          C(2*dim+i,(j-1)*dim^2+(i-1)*dim+1:(j-1)*dim^2+i*dim)=1; 
       end
    end
    
    nullC = null(C);
    
    v = ones(dim^3,1)/dim^2;
    
    new_A = A*nullC;
    new_b = zeros(dim^3,1)-A*v;
body.A = new_A;
body.b = new_b;
    body.p = zeros(size(new_A,2),1);
    volume = -1;
elseif strcmp(shape,'zonotope')==1
    %set the zonotope as the unit cube
    new_vecs = randn(dim,dim);
    for i=1:dim
        new_vecs(i,:)=new_vecs(i,:)/norm(new_vecs(i,:));
    end
    body = [eye(dim,dim); new_vecs];
    p = 0.5*ones(dim,1);
    volume = -1;
elseif strcmp(shape,'cube')==1
    body.A = [eye(dim,dim); -eye(dim,dim)];
    body.b = ones(2*dim,1);
    body.A_eq = [];
    body.b_eq = [];
    body.p = zeros(dim,1);
    volume = 2^dim;
elseif strcmp(shape,'long_box')==1
    body.A = [eye(dim,dim); -eye(dim,dim)];
    body.b = ones(2*dim,1);
    body.b(1)=100;
    body.b(dim+1)=100;
    body.A_eq = [];
    body.b_eq = [];
    body.p = zeros(dim,1);
    volume = 100*2^dim;
elseif strcmp(shape,'standard_simplex')==1
    body.A = [-eye(dim); ones(1,dim)];
    body.b = ones(dim+1,1)/(dim+1);
    body.p = zeros(dim,1);   
    volume=1/factorial(dim);
elseif strcmp(shape,'isotropic_simplex')==1
    %we need to do some computation to determine this one
    x=zeros(dim,dim+1);
    for i=1:dim
        x(i,i)=sqrt(1-norm(x(:,i))^2);
        for j=i+1:dim+1
            x(i,j)=(-1/dim-dot(x(:,i),x(:,j)))/x(i,i);
        end
    end
    
    vol_mat = zeros(dim,dim);
    for i=2:dim+1
        vol_mat(i-1,:)=x(:,i)-x(:,1);
    end
    
    body.A = x';
    body.b = ones(dim+1,1);
    body.p = zeros(dim,1);   
    volume=dim^dim/factorial(dim)*abs(det(vol_mat));
else
    fprintf('Possible shapes:\n\n');
    fprintf('cube\n');
    fprintf('standard_simplex\n');
    fprintf('isotropic_simplex\n');
    fprintf('long_box\n');
    fprintf('ball\n');
    fprintf('ellipsoid\n');
    fprintf('birkhoff\n');
    fprintf('triply_stochastic\n');
    fprintf('\nExample usage: P = makeBody(''cube'',10);\n');
    return;
end
end