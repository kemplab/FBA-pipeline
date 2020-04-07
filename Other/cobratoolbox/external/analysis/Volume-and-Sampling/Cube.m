function [volume,actual_vol] = Cube(dim, eps, type)
%CUBE This function will call Volume with a cube with the specified
%dimension/accuracy. If type=2, then the cube will be randomly
%linearly transformed. The volume of the cube is 2^dim*det(T), where T is
%the linear transformation.

[P] = makeBody('cube', dim);

if type==2
    T = randn(dim,dim);
else
    T = eye(dim);
% [T,~]=qr(randn(dim,dim));
end

P.A = P.A*T;

if type==2
    [volume] = Volume(P,[],eps,'-round 100000');
elseif type==3
    [volume] = Volume(P,[],eps,'-walk har -verb 2 -plot');
else
    [volume] = Volume(P,[],eps);
end

actual_vol = 2^dim/abs(det(T));
fprintf('Computed Volume: %e\n', volume);
fprintf('Actual Volume: %e\n', actual_vol);
fprintf('Error: %f\n', abs(actual_vol-volume)/actual_vol);

end