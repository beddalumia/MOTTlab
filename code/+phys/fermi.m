function f = fermi(w,beta)
%% FERMI-DIRAC distribution
% w is an array of real values: energies
% beta is a positive float: inverse temperature

if(beta==inf)
   fh = heaviside(-w);
end

exponential = exp(beta*w);
f = 1./(exponential+1);

if fh ~= f
   error("Wrong Heaviside!");
end

end

