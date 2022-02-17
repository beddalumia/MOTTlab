function f = fermi(w,beta)
%% FERMI-DIRAC distribution
%
% w is an array of real values: energies
% beta is a positive double: inverse temperature
%
% NB. you can totally pass 'inf' to beta and get correct T=0 results :)

    switch beta   
        
        case inf    
            
            f = heaviside(-w);  
            
        otherwise 
            
            exponential = exp(beta*w);
            f = 1./(exponential+1);   
            
    end

end

function y = heaviside(x)
         y = (x > 0) + 0.5*(x == 0);
end
