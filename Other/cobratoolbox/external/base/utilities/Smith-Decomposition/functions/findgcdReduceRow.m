function [LAC S] = findgcdReduceRow(LA,it,S)
% Work on row it + 1
gcd_found = 0;
itt = 0;
while (~gcd_found)
% for ii=1:1
    Ait = LA(it+1:end,it+1:end);
    A1 = Ait(1,:);A1NaN = A1;A1NaN(~A1NaN)=NaN;
    [min_val min_id] = min(abs(A1NaN));
    %Check if gcd is already present
    A1_eval = A1;
    A1_eval(min_id) = [];
    [g_eval c_eval d_eval] = gcd(A1_eval,A1(min_id));
    A1_fact = zeros(1,length(A1));
    for i=1:length(A1)
        facts = factor(abs(A1(i)));
        A1_fact(1,i) = facts(1);     
    end
    if sum(abs(c_eval)) == 0 && length(unique(A1_fact))~=1
        %Check if it is not reducible
            gcd_found = 1;
            break; 
    end
    for i=1:size(A1,2)
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
                if abs(qc*A1(min_id))>abs(A1(i))  
                    qc_new = abs(qc) - 1;
                    qc = sign*qc_new;
                end
                
                S = S*T3(size(LA,2),i+it,min_id+it,qc)';
                LA = LA*T3(size(LA,2),i+it,min_id+it,qc)';
                
            else
                c = -A1(i)/A1(min_id);                
                S = S*T3(size(LA,2),i+it,min_id+it,c)';
                LA = LA*T3(size(LA,2),i+it,min_id+it,c)';
                
            end
            
            Ait1 = LA(it+1:end,it+1:end);
            A11 = Ait1(:,1);A1NaN1 = A11;A1NaN1(~A1NaN1)=NaN;
            [min_val_1 min_id_1] = min(abs(A1NaN1));
            if min_val_1 == 1
                break;
            end
        end
    end
    itt=itt+1;
end

LAC = LA;
%GCD found now swap to highest diagonal position
S = S*T2(size(LA,2),it+1,min_id+it)';
LAC = LAC*T2(size(LA,2),it+1,min_id+it)';

% Set remaining entries to 0
for i=it+2:size(LA,2)
        c = LAC(it+1,i)/LAC(it+1,it+1);
        S = S*T3(size(LA,2),i,it+1,-c)';
        LAC = LAC*T3(size(LA,2),i,it+1,-c)';
end

end