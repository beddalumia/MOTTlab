function [v,Giv] = matsubara(w,Gw,beta,vcut)
%% Computes the interacting thermal Green's function
%
%           +inf
%           ⌠
%           ⎮       A(ω)    
%   G(iν) = ⎮ dω ──────────,  ν = (2n+1)πT
%           ⎮     iν  -  ω  
%           ⌡
%        -inf

    if beta > 1e4
       T = 1/1e04;
    else
       T = 1/beta;
    end

    % Define the interacting spectral function A(w)
    Aw = -imag(Gw)/pi;
    % Define the Fermionic thermal frequencies ν
    v = (pi*T):(2*pi*T):(vcut); N = length(v);
    % Compute G(iν) by Hilbert transformation
    Giv = zeros(N,1); dw = abs(w(2)-w(1));
    for n = 1:N
        dGiv = dw * Aw ./ (1i*v(n) - w);
        Giv(n)  = sum(dGiv);
    end
    
end
