function [vF,Dfig] = dirac(w,gloc)
%% Dirac velocity from an appropriate fit of the spectral function
%
%  Input:
%       w    : real valued array, ω domain
%       gloc : complex valued array, G(ω) function
%  Output:
%       vF   : sqrt( 2π/A'(ω=0+) ), A(ω) = -1/π Im G(ω)
%
%% BSD 3-Clause License
% 
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.

                                                               global DEBUG
    % Only empty states
    dos = -1/pi*imag(gloc(w>0));
    w = w(w>0);

    % Fitting domain
    w_th = 0.01; % threshold (to be tuned)
    xToFit = w(w<=w_th);

    % Windowes DOS
    yToFit = dos(w<=w_th);
    LinearModel = polyfit(xToFit,yToFit,1);

                                                                   if DEBUG
    % debug / fine-tuning of w_th
    yModel = LinearModel(1)*w + LinearModel(2);                             
    Dfig = figure("Name",'Debug on v_F determination','Visible','on');               
    plot(w,dos,':'); hold on                                         
    plot(w, yModel,'r');
    xlabel('\omega');
    legend('A(\omega)',"A'(0)\times\omega")
    xlim([0,1]);
    a = unique(min(dos));
    b = unique(max(dos));
    if a<b, ylim([a,b]); end                                
                                                                   end

    % In the proximity of the Dirac cone: A(ω)~2π/vF^2
    vF = real(sqrt(2*pi/LinearModel(1))); 

end


