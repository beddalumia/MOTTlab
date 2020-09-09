function H = semiHilbert(zeta,D)
%% HILBERT Transform, assuming a semicircular density of states of radius D
%  TO BE COMMENTED BETTER...
s = sqrt(zeta.^2 - D^2);
s = sign(imag(s)) .* s;
H = 2 * (zeta - s) / D^2;
end

