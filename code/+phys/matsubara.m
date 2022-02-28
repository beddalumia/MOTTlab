function [iw,Giw,Gtau] = matsubara(Gw,Sw,D,dos,beta)
%% Computes the noninteracting Thermal Green's function
%
%           +inf
%           ⌠
%           ⎮  dw      DOS(w)
%   G(iw) = ⎮ ────  ───────────,  iw = (2n+1)πT
%           ⎮  2π    iw  -  w
%           ⌡
%        -inf

    if isinf(beta)
       error('No zero temperature overloading for matsubara frequencies defined yet.')
    end

    % Define the Fermionic thermal frequencies
    iw = pi/beta .* (1:2:1024);
    % Compute G(iw) by Hilbert transformation
    Giw = phys.gloc(1i*iw,D,dos) / pi*2;
    % Compute G(τ) by Fourier transformation 
    Gtau = fft(Giw);

end
