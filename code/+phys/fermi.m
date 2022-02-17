function f = fermi(w,beta)
%% FERMI-DIRAC distribution
%
% w is an array of real values: energies
% beta is a positive double: inverse temperature
%
% NB. you can totally pass 'inf' to beta and get correct T=0 results :)

                                                               global DEBUG
exponential = exp(beta*w);
f = 1./(exponential+1);
                                                                   if DEBUG
if(beta==inf)
   fh = heaviside(-w);
   if any(fh ~= f)
      error("Fermi function at T=0 does not give a heaviside step!");
   end
end
                                                                   end
end


function y = heaviside(x)
         y = (x > 0) + 0.5*(x == 0);
end
