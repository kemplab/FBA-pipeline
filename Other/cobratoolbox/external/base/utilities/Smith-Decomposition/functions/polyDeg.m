function [deg] = polyDeg(poly_sym)
    poly = sym2poly(poly_sym);
    if length(poly)==1 && poly==0
        deg = inf;
    else
        deg = length(poly)-1;
    end
end