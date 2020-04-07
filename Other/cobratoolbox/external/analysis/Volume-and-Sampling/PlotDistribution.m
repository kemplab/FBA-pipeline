function [] = PlotDistribution( pts , plotHandle)
%PLOTDISTRIBUTION Plot the distribution of pts

%     np = min([its(i) P]);
    np = size(pts,1);
    h = 15;
    X = linspace(min(pts(1:np,1)),max(pts(1:np,1)),h);
    Y = linspace(min(pts(1:np,2)),max(pts(1:np,2)),h);
    Z = zeros(h,h);
    for ij=1:np
        [~,x_ind] = min(abs(X-pts(ij,1)));
        [~,y_ind] = min(abs(Y-pts(ij,2)));
        mult = 1;
        if x_ind==1 || x_ind==h
            mult = mult*2;
        end
        if y_ind==1 || y_ind ==h
            mult = mult*2;
        end
        Z(y_ind,x_ind) = Z(y_ind,x_ind)+mult/np*(h-1)^2;
    end
    
    figure(plotHandle);
    surf(X,Y,Z);
    drawnow;
    axis([X(1) X(end) Y(1) Y(end) 0 1.25*max(max(Z))]);

end