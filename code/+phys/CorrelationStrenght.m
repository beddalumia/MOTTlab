%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

function S = CorrelationStrenght(w,sloc)
%% Computes the Strenght of Correlations, as defined in PRL 114 185701
%  Input:
%       w           : real valued array, \omega domain
%       sloc        : complex valued array, \Sigma(\omega) function
    infty = length(w);
    zero  = round(infty/2)+1;
    S = norm(sloc(zero)-sloc(infty));
end

