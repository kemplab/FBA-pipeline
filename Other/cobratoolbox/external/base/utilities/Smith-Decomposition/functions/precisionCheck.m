function [A] = precisionCheck(A)
%precision
symv = symvar(A);

    for ii=1:size(A,1)
        for jj=1:size(A,2)
                poly_coeffs = sym2poly(A(ii,jj));
                for kk=1:length(poly_coeffs)
                    if abs(poly_coeffs(kk)) < 10^-10
                        poly_coeffs(kk) = 0;
                    end
                end
                if ~isempty(symv)
                    A(ii,jj) = poly2sym(poly_coeffs,symv);
                else
                    A(ii,jj) = poly2sym(poly_coeffs);
                end
                    
        end
    end
end