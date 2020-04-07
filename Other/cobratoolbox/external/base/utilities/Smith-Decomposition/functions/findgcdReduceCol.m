function [T LA] = findgcdReduceCol(A,it,T)
%Apply row operations to reduce the column until gcd of all entries is
%found

% Work on column it + 1
gcd_found = 0;
LA = A;
itt = 0;
while (~gcd_found)
% for ii=1:5
    Ait = LA(it+1:end,it+1:end);
    A1 = Ait(:,1);A1NaN = A1;A1NaN(~A1NaN)=NaN;
    [min_val min_id] = min(abs(A1NaN));
    
    %Check if gcd is already present
    A1_eval = A1;
    A1_eval(min_id) = [];
    [g_eval c_eval d_eval] = gcd(A1_eval,A1(min_id));
    if sum(abs(c_eval)) == 0
        gcd_found = 1;
        break;
    end
    
    for i=1:size(A1,1)
            if i~=min_id && A1(i)~=0
                [g c d] = gcd(A1(i),A1(min_id));
                if c~=0
                        q = d/c;
                        if q < 0
                            sign = -1;
                        else
                            sign = 1;
                        end
                        if abs(floor(abs(q))-abs(q))>0.4 && itt==0
                            qc = sign*(floor(abs(q)) + 1);
                        else
                            qc = sign*(floor(abs(q)));
                        end
                        
                        %Check factor
                        if abs(qc*A1(min_id))>abs(A1(i))  && itt>0
                            qc_new = abs(qc) - 1;
                            qc = sign*qc_new;
                        end                       
                        T = T3(size(A,1),i+it,min_id+it,qc)*T;                     
                        LA = T3(size(A,1),i+it,min_id+it,qc)*LA;                                           
                else
                    c = A1(i)/A1(min_id);
                    if c>0
                        c=-c;
                    end
                    T = T3(size(A,1),i+it,min_id+it,c)*T;                 
                    LA = T3(size(A,1),i+it,min_id+it,c)*LA;                    
                end                
                  
            end
                Ait1 = LA(it+1:end,it+1:end);
                A11 = Ait1(:,1);A1NaN1 = A11;A1NaN1(~A1NaN1)=NaN;
                [min_val_1 min_id_1] = min(abs(A1NaN1));
                if min_val_1 == 1 %lowest order gcd or lowest poly order
                    break;
                end 
    end
    itt=itt+1;
end

%GCD found now swap to highest diagonal position
if min_id~=1
    T = T2(size(A,1),it+1,min_id+it)*T;
    LA = T2(size(A,1),it+1,min_id+it)*LA;
end
% Change sign of gcd if negative
if LA(it+1,it+1)<0
    T = T1(size(A,1),it+1,-1)*T;
   LA = T1(size(A,1),it+1,-1)*LA;
end

% Set remaining entries to 0
for i=it+2:size(A,1)
        c = LA(i,it+1)/LA(it+1,it+1);
        T = T3(size(A,1),i,it+1,-c)*T;
        LA = T3(size(A,1),i,it+1,-c)*LA;
end

end