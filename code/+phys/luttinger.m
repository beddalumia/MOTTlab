function [I,Lfig] = luttinger(w,sloc,gloc,axis)
%% Computes the Luttinger integral, either on the real or imaginary axis
%
%  Input:
%       w     : real valued array, ω domain OR iω matsubara discretization
%       sloc  : complex valued array, ∑(ω) OR ∑(iω) function
%       gloc  : complex valued array, G(ω) OR G(iω) function
%       axis  : optional string, enforcing algorithm ['w' OR 'iw']
%  Output:
%       I = \frac{2}{\pi}\Im\int_{-\infty}^{0} dw G_loc(w) d\Sigma(w) / dw
%         = \frac{1}{\pi}\Im\int_{-\infty}^{\infty}dwG_loc(w)d\Sigma(w)/dw
%
%% BSD 3-Clause License
% 
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.
                                                               global DEBUG
    if nargin < 4     % Automatic detection of the axis,
       if all(w>0)    % based on reasonable assumptions.
          axis = 'matsubara';
       else
          axis = 'real-axis';
       end
    end
    
    switch axis
        
        case {'real','real-axis','spectral','w'}
    
            % w = w(w<=0);          % Faster (half points to sum)
            % sloc = sloc(w<=0);    % but apparently less accurate
            % gloc = gloc(w<=0);    % --> better to exploit ph-sym

            ds  = diff(sloc); % dSigma -> already includes dw!
            g_  = gloc(1:end-1);
            integrand  = imag(g_.*ds);
            I = 1/pi*sum(integrand);
            I = abs(I); % J. Phys.: Condens. Matter 28 (2016) 025601
if DEBUG
            Lfig = figure("Name",'Luttinger integrand','Visible','off');
            plot(w(1:end-1),integrand);
            xlabel('\omega');
            ylabel('Im[G(\omega)\partial\Sigma/\partial\omega]');
            ylim([-0.1,0.2]);
end                                                                  
                                                                   
        case {'imag','imag-axis','thermal','matsubara','iw','iv'}
            
            ds  = diff(sloc); % dSigma -> already includes dw!
            g_  = gloc(1:end-1);
            integrand  = real(g_.*ds);
            I = 2*w(1)/pi*sum(integrand);
            I = abs(I); % J. Phys.: Condens. Matter 28 (2016) 025601
if DEBUG
            Lfig = figure("Name",'Luttinger integrand','Visible','off');
            plot(w,imag(gloc)); hold on; plot(w,imag(sloc));
            plot(w(1:end-1),integrand);
            xlabel('i\omega');
            ylabel('Re[G(i\omega)\partial\Sigma/\partial i\omega]');
end  
           
        otherwise
            
    end
    
end


