%where the points are weighted by the distribution exp(-a_i||x||^2)
function [new_pt] = rand_exp_range_coord(l,u,a_i)
sigma = sqrt(1/2/a_i);

if l > 3*sigma || u < -3*sigma
    if l > 3*sigma
        if (u-l) <= 2*sigma
           %do simple rejection sampling
           while 1
              new_pt = rand()*(u-l) + l;
              h = rand();
              if h <= exp((l^2-new_pt^2)/(2*sigma^2))
                 break; 
              end
           end
        else
            %otherwise sample from exponential
            mu = (2*sigma^2)/l;

            while 1
               new_pt = exprnd(mu) + l;
               
%                l
%                u
%                sigma
%                cutoff = exp(-new_pt^2/(2*sigma^2) + new_pt/mu)
               if new_pt <= u && rand() <= exp(-new_pt^2/(2*sigma^2) + new_pt/mu)
                   break;
               end
            end
        end
    else
        if (u-l) <= 2*sigma
            while 1
                new_pt = rand()*(u-l)+l;
                h = rand();
                if h<=exp((u^2-new_pt^2)/(2*sigma^2))
                    break;
                end
            end
            
        else
            %otherwise sample from exponential
            mu = -2*sigma^2/u;
            while 1
               new_pt = u-exprnd(mu);
               if new_pt >= l && rand() <= exp(-new_pt^2/(2*sigma^2)-new_pt/mu)
                   break;
               end
            end
        end
    end

elseif a_i>1e-8 && u-l>=2/sqrt(2*a_i)
    %select from the 1d Gaussian chord if enough weight will be inside
    %K
    
    while 1
        %sample from Gaussian along chord, and accept if inside (u,v)
        
        rn = randn(1)/sqrt(2*a_i);
        if rn>=l && rn <=u
            break;
        end
    end
    new_pt = rn;
else
    %otherwise do simple rejection sampling by a bounding rectangle
    M = get_max_coord(l,u,a_i);
    done = 0;
    %     its = 0;
    while ~done
        %         its = its+1;
        rn = rand();
        pt = (1-rn)*l+rn*u;
        r_val = M*rand();
        fn = eval_exp(pt, a_i);
        if r_val<fn
            done = 1;
            new_pt = pt;
        end
    end
end
end

function [ret1] = get_max_coord(l, u, a_i)
%get the maximum value along the chord, which is the height of the bounding
%box for rejection sampling
if l<0 && u > 0
    ret1 = 1;
else
    ret1 = max(eval_exp(l,a_i),eval_exp(u,a_i));
end

end
