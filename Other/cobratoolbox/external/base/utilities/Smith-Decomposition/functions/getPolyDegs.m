function[polydegs] = getPolyDegs(Ai)
    polydegs = zeros(length(Ai),1);
    for ii=1:length(Ai)
        polydegs(ii,1) = polyDeg(Ai(ii));
    end
end