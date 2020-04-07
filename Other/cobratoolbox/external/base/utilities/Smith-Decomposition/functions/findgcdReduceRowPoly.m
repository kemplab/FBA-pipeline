function [LAC S] = findgcdReduceRowPoly(LA,it,S)
%Apply column operations (S1=T1',S2=T2',S3=T3') to reduce the row until 
%gcd of all entries is found

gcd_found = 0;
LAC = LA;
symv = symvar(LA);
Ait = LAC(it+1:end,it+1:end);
    
if ~(isempty(Ait) || isempty(symvar(Ait)))
    while (~gcd_found)
        Ait = LAC(it+1:end,it+1:end);
        A1 = Ait(1,:);
        A1 = precisionCheck(A1);
        LAC(it+1:end,it+1:end) = precisionCheck(LAC(it+1:end,it+1:end));
        
        [min_deg min_id] = minPolyDeg(A1);
        polydegs = getPolyDegs(A1);      
        
        if min_deg == 0
            gcd_found = 1;
            break;
        else
            [qs rs] = polyDivs(A1,Ait(1,min_id));
        end
        
        if sum(abs(rs))==0
            gcd_found = 1;
            break;
        end
        for i = 1:size(A1,2)
            if i~=min_id && A1(i)~=0
            [q,r] = deconv(sym2poly(A1(i)),sym2poly(A1(min_id)));
                if length(q)==1 %quotient is constant
                    alpha = -q(1);
                    S = S*T3(size(LA,2),min_id+it,i+it,alpha);
                    LAC = LAC*T3(size(LA,2),min_id+it,i+it,alpha);    
                else %quotient is a polynomial
                    alpha = -poly2sym(q,symv);
                    S = S*T3(size(LA,2),min_id+it,i+it,alpha);
                    LAC = LAC*T3(size(LA,2),min_id+it,i+it,alpha);    
                end                       
            end
        end
    end 

    % If GCD not in highest diagonal position --> swap!
    if min_id~=1
%         display('swap to highest')
        S = S*T2(size(LA,2),min_id+it,it+1);
        LAC = LAC*T2(size(LA,2),min_id+it,it+1);
    end

    gcd_poly = sym2poly(LAC(it+1,it+1));
   
    % If GCD is a non-monic polynomial  -->  make monic!
    if length(gcd_poly)>1 && gcd_poly(1)~=1 
%         display('make monic')
        S = S*T1(size(LA,2),it+1,1/gcd_poly(1));
        LAC = LAC*T1(size(LA,2),it+1,1/gcd_poly(1));
    end
    
    % If GCD is a non-monic  -->  make monic!
    if length(gcd_poly)==1 && gcd_poly(1)~=1
%         display('reduce to 1')
        S = S*T1(size(LA,2),it+1,1/gcd_poly(1));
        LAC = LAC*T1(size(LA,2),it+1,1/gcd_poly(1));
        %precision
        LAC(it+1,it+1) = double(LAC(it+1,it+1));
    end

    % If GCD negative --> negate!
    gcd_poly = sym2poly(LAC(it+1,it+1));
    if gcd_poly(1)<0
%         display('negate')
        T1(size(LA,2),it+1,-1)
        S = S*T1(size(LA,2),it+1,-1);
        LAC = LAC*T1(size(LA,2),it+1,-1);
    end
    
    % Set remaining entries to 0
    for i=it+2:size(LA,2)
%         display('remaining 0s')
        if ~(length(LAC(it+1,i))==1 && LAC(it+1,i)==0)
            [q r] = deconv(sym2poly(LAC(it+1,i)),sym2poly(LAC(it+1,it+1)));
            if length(q)==1 %quotient is constant
                c = -q(1);
                if c == 0
                    c = -r;
                end
            else %quotient is a polynomial
                c = -poly2sym(q,symv);
            end
            S = S*T3(size(LA,2),it+1,i,c);
            LAC = LAC*T3(size(LA,2),it+1,i,c);
            %precision
            poly = sym2poly(LAC(it+1,i));
            for kk=1:length(poly)
                if abs(poly(kk)) < 10^-10
                    poly(kk) = 0;
                end
            end
            LAC(it+1,i) = poly2sym(poly,symv);
        end
    end
end

end