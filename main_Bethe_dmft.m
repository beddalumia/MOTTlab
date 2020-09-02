%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clear all; clc;
%% INPUT: Physical Parameters 
D = 1;                  % Bandwidth
U = 2.7;                % On-site Repulsion    } Overriden if PhaseDiagram
beta = 50;             % Inverse Temperature  } flag is set to true...

%% INPUT: Boolean Flags
InsulatingBIAS = 1;     % Changes initial guess of gloc (should favour Mott phase)
PhaseDiagram = 1;       % Overrides given U and beta values and performs a full scan
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
w = linspace(-6, 6, Resolution); 

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
    i = 0; T = 10^(-6);
    while T < 0.05
        i = i + 1;
        j = 0; 
        U = 0; 
        while U <= 5.00  
            j = j + 1; beta = 1/T;
            [gloc{i,j},Sigma_loc{i,j}] = DMFT_loop(gloc_0,w,D,U,beta,loops,mixing,error);
            U = U + Ustep;
        end
        T = T + Tstep; 
    end
    
else
    %% Single (U, beta) pair calculation
    [gloc,Sigma_loc] = DMFT_loop(gloc_0,w,D,U,beta,loops,mixing,error); 
end

%% Plotting stuff (calling a script...)
plotter_Bethe_dmft;
