function c = fconv(a, b, shape)
%% Fast convolution of long arrays, exploiting NEXTPOW2 and built-in FFTW wrappers
%
%   C = FCONV(A, B, SHAPE) convolves vectors A and B. The length of the resulting
%                          vector is determined by the optional SHAPE tag.
%
%     'full'  - (default) returns the full convolution: 
%               LENGTH(C) = MAX([LENGTH(A)+LENGTH(B)-1,LENGTH(A),LENGTH(B)])
%     'same'  - returns the central part of the convolution:
%               LENGTH(C) == LENGTH(A)
%     'valid' - returns only those parts of the convolution 
%               that are computed without the zero-padded edges: 
%               LENGTH(C) is MAX(LENGTH(A)-MAX(0,LENGTH(B)-1),0)
%
%   The FCONV implementation is much faster(*) than standard generic CONV if:
%
%       - A and B are very large
%       - A and B have similar length
%       - SHAPE is different from 'valid' (x)
%
%       > if length(A) ~= length(B) by a good amount, you should call indeed
%         the standard CONV function, for it may be much faster(-) than FCONV.
% 
%   See also CONV, NEXTPOW2, FFT, IFFT.
%
%% NOTES:
%
%  (*) I have seen speedups up to 40x for length(A)=length(B)=2^16.
%
%  (x) if shape == 'valid' this function just wraps the built-in CONV
%
%  (-) I have seen slowdowns up to 50x if lenght(A)-lenght(B)
%      is a very large multiple of the smaller between the two.
%
%
%% BSD 3-Clause License
%
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.

 if nargin < 3
    shape = 'full';                      % default shape
 elseif isstring(shape)
    shape = char(shape);
 end

 if shape(1) == 'v' || shape(1) == 'V'   %  shape 'valid'
    c = conv(a, b, shape); return
 end
 
 ABeq = false;
 if a == b
     ABeq = true;
 end
 
 ABpos = false;
 if isreal(a) && isreal(b)
    if all(a>=0) && all(b>=0)
       ABpos = true; 
    end
 end
    
    Na = length(a);
    Nb = length(b);
    
    Nc = Na + Nb - 1;           % Standard full convolution length
    Npow = pow2(nextpow2(Nc));  % FFT OPTIMIZATION TRICK: see doc nextpow2
    a = fft(a, Npow);           % Fastest Fourier transform in the West
 if ABeq
    b = a;                      % Avoid second FFT if init a == init b
 else
    b = fft(b, Npow);           % Compute second FFT otherwise
 end
    c = a .* b;                 % Fast vectorized scalar product
    c = real(ifft(c, Npow));    % FFTW's inverse fast Fourier transform
    
 if ABpos       
    c = abs(c);                 % Enforce semi-definite positive property
 end
    
 if shape(1) == 'f' || shape(1) == 'F'   % shape 'full'
    c = c(1:Nc); return                    
 end   
    
 if shape(1) == 's' || shape(1) == 'S'   % shape 'same'
    c = c(1:Nc);
    Ncut = (Nc-Na)/2;
    c = c(1+ceil(Ncut):end-floor(Ncut)); return                        
 end
   
end



