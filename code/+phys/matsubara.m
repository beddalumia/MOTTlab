function [v,Fiv] = matsubara(w,Fw,beta,vcut,bcut)
%% MATSUBARA representation from spectral decomposition.
%
%  Input:
%       w           : real valued array, ω domain, must be centrosymmetric
%       Fw          : complex valued array, F(ω) function, must be causal
%       beta        : real double, inverse temperature, must be positive
%       vcut        : real double, cutoff on iν-axis, must be positive
%       bcut        : real double, cutoff on beta [OPTIONAL, default: 1e4]
%  Output:
%       v           : real valued array, positive matsubara frequencies ν
%       Fiv         : complex valued array, F(iν) function, defined as:
%
%                               +inf
%                               ⌠                 
%                               ⎮       A(ω)      ⎧ A(ω) = -Im[F(ω)]/π
%                       F(iν) = ⎮ dω ──────────,  ⎨
%                               ⎮     iν  -  ω    ⎩ ν = (2n+1)πT
%                               ⌡                   
%                            -inf
%
%   NB. F(•) would typically be either Gloc(•) or Sloc(•)
%
%% BSD 3-Clause License
%
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.

    if nargin < 5
       bcut = 1e3; 
    end

    if beta > bcut
       fprintf(2,'WARNING:\n');
       fprintf(2,'Analytic continuation to thermal axis\n');
       fprintf(2,'performed with fictitious beta = %d\n\n',bcut);
       T = 1/bcut;
    else
       T = 1/beta;
    end
    
    % Define the interacting spectral function A(ω)
    Aw = -imag(Fw)/pi;
    
    % Define the Fermionic thermal frequencies ν
    if pi*T < vcut
       v = (pi*T):(2*pi*T):(vcut); N = length(v);
    else
       error('Whooops, vcut < πT, matsubara aborting!'); 
    end
    
    % Compute F(iν) by Hilbert transformation
    
    %% DIRECT SUM-OVER-ALL-W IN A IW LOOP
    Fiv = zeros(1,N); dw = abs(w(2)-w(1));
    for n = 1:N
        dFiv = dw * Aw ./ (1i*v(n) - w);
        Fiv(n) = sum(dFiv);
    end
    
    Ftest = Fiv;
    
    %% HANDLE-VECTORIZED SUM-OVER-ALL-W
    dw = abs(w(2)-w(1));
    dFiv = @(z) dw * Aw' ./ (z - w');
    Fiv = sum(dFiv(1i*v));
    
    err = norm(Fiv-Ftest)/norm(Ftest);
    if err > 10*eps
       warning('problems with vectorization? err = %.16f',err); 
    end
    
end



