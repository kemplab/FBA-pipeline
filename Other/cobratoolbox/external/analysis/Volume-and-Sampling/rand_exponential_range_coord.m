%where the points are weighted by the distribution exp(-a_i||x||)
function [new_pt] = rand_exponential_range_coord(l,u,a)


ratio = exp(-a*l)/exp(-a*u);
ratio = max(ratio,1/ratio);

if ratio < 4
   %do rejection sampling
   max_val = max(exp(-a*l), exp(-a*u));
   while 1
      new_pt = rand()*(u-l) + l;
      if rand()*max_val <= exp(-a*new_pt)
          break;
      end
   end
else
    if a<0
        flipped = 1;
        a = -a;
    else
        flipped = 0;
    end

    %sample from the exponential
%     if exp(-a*l)<exp(-a*u)
%         flipped = 1;
%     else
%         flipped = 0;
%     end
    tries = 0;
    while 1
        tries = tries+1;
%        r = exprnd(1/a);
       new_pt = l + exprnd(1/a);
       if l <= new_pt && new_pt <= u
           break;
       end
       if tries > 200
          fprintf('uh-oh'); 
       end
       if flipped
           new_pt = u-(new_pt-l);
       end
    end
end
% 
% if a_i>1e-8 && u-l>=2/sqrt(2*a_i)
%     %select from the 1d Gaussian chord if enough weight will be inside
%     %K
%     
%     while 1
%         %sample from Gaussian along chord, and accept if inside (u,v)
%         
%         rn = randn(1)/sqrt(2*a_i);
%         if rn>=l && rn <=u
%             break;
%         end
%     end
%     new_pt = rn;
% else
%     %otherwise do simple rejection sampling by a bounding rectangle
%     M = get_max_coord(l,u,a_i);
%     done = 0;
%     %     its = 0;
%     while ~done
%         %         its = its+1;
%         rn = rand();
%         pt = (1-rn)*l+rn*u;
%         r_val = M*rand();
%         fn = eval_exp(pt, a_i);
%         if r_val<fn
%             done = 1;
%             new_pt = pt;
%         end
%     end
% end
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
