function [T LA] = findgcdReduceColPoly(A,it,T)
%Apply row operations (T1,T2,T3) to reduce the column until 
%gcd of all entries is found

gcd_found = 0;
LA = A;
symv = symvar(A);
Ait = LA(it+1:end,it+1:end);
A1 = Ait(:,1);

if ~(isempty(Ait) || (sum(abs(A1))==0))
    while (~gcd_found)
        Ait = LA(it+1:end,it+1:end);
        A1 = Ait(:,1);            
        A1 = precisionCheck(A1);
        LA(it+1:end,it+1:end) = precisionCheck(LA(it+1:end,it+1:end));

        [min_deg min_id] = minPolyDeg(A1);
        polydegs = getPolyDegs(A1);
        if min_deg == 0
            gcd_found = 1;
            break;
        else
            [qs rs] = polyDivs(A1,Ait(min_id,1));
        end
        if sum(abs(rs))==0
            gcd_found = 1;
            break;
        end

        for i = 1:size(A1,1)
            if i~=min_id && A1(i)~=0
            [q,r] = deconv(sym2poly(A1(i)),sym2poly(A1(min_id)));
                if length(q)==1 %quotient is constant
                    alpha = -q(1);
                    T = T3(size(A,1),i+it,min_id+it,alpha)*T;
                    LA = T3(size(A,1),i+it,min_id+it,alpha)*LA;
                else %quotient is a polynomial
                    alpha = -poly2sym(q,symv);
                    T = T3(size(A,1),i+it,min_id+it,alpha)*T;
                    LA = T3(size(A,1),i+it,min_id+it,alpha)*LA;
                end
            end
        end
    end

    % If GCD not in highest diagonal position --> swap!
    if min_id~=1
%         display('swap to highest')
        T = T2(size(A,1),it+1,min_id+it)*T;
        LA = T2(size(A,1),it+1,min_id+it)*LA;
    end
    
    % Checking precision errors
    gcd_poly = sym2poly(LA(it+1,it+1));
    
    % If GCD is a non-monic polynomial  -->  make monic!
    if length(gcd_poly)>1 && gcd_poly(1)~=1 && length(A1)~=1
%         display('make monic')
        T = T1(size(A,1),it+1,1/gcd_poly(1))*T;
        LA = T1(size(A,1),it+1,1/gcd_poly(1))*LA;
    end

    % % If GCD is a constant ~=1 --> reduce!
    if length(gcd_poly)==1 && gcd_poly(1)~=1
%         display('reduce to 1')
        T = T1(size(A,1),it+1,1/gcd_poly(1))*T;
        LA = T1(size(A,1),it+1,1/gcd_poly(1))*LA;
        %precision
        LA(it+1,it+1) = double(LA(it+1,it+1));
    end

    % Set remaining entries to 0
    for i=it+2:size(A,1)
%             display('remaining 0s')
            [q r] = deconv(sym2poly(LA(i,it+1)),sym2poly(LA(it+1,it+1)));
            if length(q)==1 %quotient is constant
                c = -q(1);
                if c == 0
                    c = -r;
                end
            else %quotient is a polynomial
                c = -poly2sym(q,symv);
            end
            T = T3(size(A,1),i,it+1,c)*T;
            LA = T3(size(A,1),i,it+1,c)*LA;
            %precision
            poly = sym2poly(LA(i,it+1));
            for kk=1:length(poly)
                if abs(poly(kk)) < 10^-10
                    poly(kk) = 0;
                end
            end
            LA(i,it+1) = poly2sym(poly,symv);

    end

end

end