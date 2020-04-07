function [min_deg min_id] = minPolyDeg(Ai)
    polydegs = zeros(length(Ai),1);
    for ii=1:length(Ai)
        polydegs(ii,1) = polyDeg(Ai(ii));
    end
    [min_deg min_id] = min(polydegs);
end