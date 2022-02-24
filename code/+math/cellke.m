function [k,e] = cellke(m,tol)
%% CELLKE Complete elliptic integral. Allowing complex input.
%   [K,E] = CELLKE(M) returns the value of the complete elliptic
%   integrals of the first and second kinds, evaluated for each
%   element of M. Domain: M ∈ ℂ (see [2] for details).
%   
%   [K,E] = CELLKE(M,TOL) computes the complete elliptic integrals to
%   the accuracy TOL instead of the default TOL = EPS(CLASS(M)).
%
%   Some definitions of the complete elliptic integrals use the modulus
%   k instead of the parameter M. They are related by M = k^2.
%
%   Class support for input M:
%      float: double, single
%
%   CELLKE extends the native implementation provided by ELLIPKE,
%   generalizing the method of the arithmetic-geometric mean [1],
%   so to allow complex values of M, as described in [2].
%
%   The algorithm has been tested against ELLIPTICK, included in 
%   the SYMBOLIC MATH TOOLBOX: the returned values appear to match
%   to machine precision, with an average ≈300x speedup.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications, 1965, 17.6.
%   [2] B.C. Carlson, "Numerical computation of real or complex 
%       elliptic integrals" Numerical Algorithms volume 10, 
%       pages 13–26 (1995) [arXiv:math/9409227]
%
%   See also ELLIPKE, ELLIPTICK.
%
%% BSD 3-Clause License
%
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.

if nargin<1
  error('Not enough input arguments.'); 
end

classin = superiorfloat(m);

if nargin<2, tol = eps(classin); end

if isempty(m), k = zeros(size(m),classin); e = k; return, end

if ~isreal(tol)
    error('Second argument TOL must be real.');
end

if ~isscalar(tol) || tol < 0 || ~isfinite(tol)
  error('Second argument TOL must be a finite nonnegative scalar.');
end

a0 = 1;
b0 = sqrt(1-m);
c0 = NaN;
s0 = m;
i1 = 0; mm = Inf;
while mm > tol
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    i1 = i1 + 1;
    w1 = 2^i1*norm(c1).^2;
    mm = max(w1(:));
    
    % Test for stagnation (may happen for TOL < machine precision)
    if isequal(c0, c1)
        error('CELLKE did not converge. Consider increasing TOL.');
    end
    
    s0 = s0 + w1;  
    a0 = a1;  b0 = b1;  c0 = c1;
end
k = pi./(2*a1);
e = k.*(1-s0/2);
im = find(m==1);
if ~isempty(im)
    e(im) = ones(length(im),1);
    k(im) = inf;
end
