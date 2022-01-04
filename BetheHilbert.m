function gloc = BetheHilbert(zeta,D)
%% HILBERT transform of the BETHE semicircular DOS of radius D
%  TO BE COMMENTED BETTER...
s    = sqrt(zeta.^2 - D^2);
s    = sign(imag(s))  .* s;
gloc = 2 * (zeta - s) / D^2;
end

