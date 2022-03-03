function Z = zetaweight(w,sloc,axis)
%% ZETAWEIGHT extracts the quasiparticle weight from either ∑(ω) or ∑(iω)
%
%  Input:
%       w    : real valued array, ω domain OR iω matsubara discretization
%       sloc : complex valued array, ∑(ω) OR ∑(iω) function
%       axis : optional string, enforcing algorithm ['w' OR 'iw']
%  Output:
%       Z    : 1/(1-d/dω[Re∑(ω=0)]) OR 1/(1-Im∑(iπT)/πT)
%
%% BSD 3-Clause License
% 
%  Copyright (c) 2020-2022, Gabriele Bellomia
%  All rights reserved.

    if nargin < 3               % Automatic detection of the axis,
       if all(w>0)              % based on reasonable assumptions.
          axis = 'matsubara';
       else
          axis = 'real-axis';
       end
    end

    switch axis
        
        case {'real','real-axis','spectral','w'}

             % Hardcoded threshold to define proper fitting window
             w_th = 1e-2; % (don't base it on wres, it's dangerous)
             % Frequency windowing of ω and ∑(ω)
             xToFit = w(abs(w)<=w_th);
             yToFit = real(sloc(abs(w)<=w_th));
             % Linear fitting the windowed data
             LinearModel = polyfit(xToFit,yToFit,1);
             % Nozieres Theorem: Re∑(ω) = A + (1-1/Z)ω + O(ω^3)
             Z = 1/(1-LinearModel(1)); 
             
        case {'imag','imag-axis','thermal','matsubara','iw','iv'}
            
             % Nozieres Theorem: Im∑(iω) = A + (1-1/Z)iω + O(iω^3)
             Z = 1/(1-imag(sloc(1))/w(1));  % Recall: ω = (2n+1)πT
            
        otherwise
            
             warning('Invalid axis: zetaweight will skip computation.');
             Z = 0;
             
    end

end
