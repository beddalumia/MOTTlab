%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

clear all; clc;
%% INPUT: Physical Parameters 
D = 1;                  % Bandwidth
U = 5;                  % On-site Repulsion    } Overriden if PhaseDiagram
beta = 50;              % Inverse Temperature  } flag is set to true...

%% INPUT: Boolean Flags
InsulatingBIAS = 0;     % Changes initial guess of gloc (should favour Mott phase)
PhaseDiagram = 0;       % Overrides given U and beta values and performs a full scan
ZeroTempSpan = 0;       % If both PhaseDiagram and ZeroTemp are set to true the scan
                        % spans only U values, with beta fixed to 1000.

%% INPUT: Technical Parameters
loops = 1000;           % Max number of DMFT iterations 
error = 10^(-5);        % Convergence threshold for self-consistency
mixing = 0.1;           % Mixing parameter for DMFT iterations (=1 means full update)
Resolution = 2^12;      % Energy resolution in real-frequency axis
Ustep = 0.05;           % Hubbard U incremental step for Phase Diagrams
Tstep = 0.005;          % Temperature incremental step for Phase Diagrams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main

% Frequency Values
w = linspace(-6, 6, Resolution*10); 

% Initial guess for the local Green's function
gloc_0 = semiHilbert(w + 10^(-3)*1i,D); % D is the DOS "radius"
if InsulatingBIAS == true
   gloc_0 = 0; % i.e. No BATH :)
end


if PhaseDiagram == true && ZeroTempSpan == true
    %% Zero Temperature MIT span
    U = 0; beta = 10^6; i = 0;
    while U <= 5.00 
        i = i + 1;
        [gloc{i},Sigma_loc{i}] = DMFT_loop(gloc_0,w,D,U,beta,loops,mixing,error);         
        U = U + Ustep;
    end

elseif PhaseDiagram == true && ZeroTempSpan == false
    %% Full Phase Diagram span
    i = 0; T = 10^(-6); restart_gloc = gloc_0;
    while T < 0.05
        i = i + 1;
        j = 0; 
        U = 0; 
        while U <= 5.00  
            j = j + 1; beta = 1/T;
            [gloc{i,j},Sigma_loc{i,j}] = DMFT_loop(restart_gloc,w,D,U,beta,loops,mixing,error);
            restart_gloc = gloc{i,j};
            U = U + Ustep;
        end
        T = T + Tstep; 
    end
    
else
    %% Single (U, beta) pair calculation
    [gloc,Sigma_loc] = DMFT_loop(gloc_0,w,D,U,beta,loops,mixing,error); 
    I = LuttingerIntegral(w,Sigma_loc,gloc)
end

%% Plotting stuff (calling a script...)
%plotter_Bethe_dmft;
