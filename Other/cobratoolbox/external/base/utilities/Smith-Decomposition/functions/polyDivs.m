function [q r] = polyDivs(Ai,b)

    for i=1:length(Ai)
            [q_poly r_poly] = deconv(sym2poly(Ai(i)),sym2poly(b));
            q(i,1) = poly2sym(q_poly);
            r(i,1) = poly2sym(r_poly);
    end

end