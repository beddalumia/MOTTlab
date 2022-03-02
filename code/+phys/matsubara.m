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
    Aw = -imag(Gw)/pi; Nreal = length(Aw);
    % Define the Fermionic thermal frequencies ν
    v = (pi*T):(2*pi*T):(vcut); Nmats = length(v);
    % Compute G(iν) by Hilbert transformation
    Giv = zeros(Nmats,1); dw = abs(w(2)-w(1));
    for n = 1:Nmats
        dGiv = dw * Aw ./ (1i*v(n) - w);
        Giv(n)  = sum(dGiv);
    end
    Giv = Giv.';
    % Vectorized form
    matV = repmat(v,Nreal,1);
    matA = repmat(Aw',1,Nmats);
    matW = repmat(w',1,Nmats);
    dw   = abs(w(2)-w(1));
    dGiv = dw * matA ./ (1i*matV - matW);
    Gvec = sum(dGiv);
    
    if any(abs(Gvec-Giv)>eps)
       warning('Vectorization reduces accuracy')
       norm(Gvec-Giv)
    end
    
end
