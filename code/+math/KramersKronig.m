function K = KramersKronig(F)
%% KramersKronig transform, exploiting NEXTPOW2 and built-in FFTW wrappers
%
%       K[F(..)](:) = KramersKronig(F), built through the chain:
%
%                   >   S[F] = hilbert(F) [Analytic-Signal]
%
%                   >   H[S] = imag(S)    [Hilbert-Transform]
%
%                   >   K[H] = -H         [Kramers-Kronig]
%
% See also NEXTPOW2, HILBERT, FFT
%
%% Theoretical Background at:
%
%  https://en.wikipedia.org/wiki/Kramers–Kronig_relations
%
%  https://en.wikipedia.org/wiki/Hilbert_transform
%
%% BSD 3-Clause License
%
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.

    N = length(F)*4;        % Nyquist condition for the relevant FFTs
    P = pow2(nextpow2(N));  % FFT OPTIMIZATION TRICK: run < doc nextpow2 >
    S = hilbert(F,P);       % Wrapper of the FFTWs:   run < open hilbert >
    K = imag(-S(1:N/4));    % H = Im(AnalyticSignal): run < help hilbert >
    
end
