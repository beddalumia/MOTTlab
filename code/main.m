%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

clearvars; clc;
global DEBUG FAST

%% INPUT: Physical Parameters 
D    = 1;               % Bandwidth
U    = 1;               % On-site Repulsion
beta = 1e6;             % Inverse Temperature

%% INPUT: Boolean Flags
MottBIAS     = 0;       % Changes initial guess of gloc (strongly favours Mott phase)
Uline        = 0;       % Takes and fixes the given beta value and performs a U-driven line
Tline        = 0;       % Takes and fixes the given U value and performs a T-driven line
UTscan       = 0;       % Ignores both given U and beta values and builds a full phase diagram
DoSPECTRAL   = 1;       % Controls plotting of spectral functions
DoPLOT       = 1;       % Controls plotting of *all static* figures
DoGIF        = 0;       % Controls plotting of *animated* figures
DEBUG        = 0;       % Activates debug prints / plots / operations
FAST         = 1;       % Activates fast FFTW-based convolutions

%% INPUT: Control Parameters
mloop = 1000;           % Max number of DMFT iterations 
err   = 1e-5;           % Convergence threshold for self-consistency
mix   = 0.10;           % Mixing parameter for DMFT iterations (=1 means full update)
wres  = 2^15;           % Energy resolution in real-frequency axis
wcut  = 6.00;           % Energy cutoff in real-frequency axis
Umin  = 0.00;           % Hubbard U minimum value for phase diagrams
Ustep = 0.09;           % Hubbard U incremental step for phase diagrams
Umax  = 6.00;           % Hubbard U maximum value for phase diagrams
Tmin  = 5e-3;           % Temperature U minimum value for phase diagrams
Tstep = 5e-3;           % Temperature incremental step for phase diagrams
Tmax  = 5e-2;           % Temperature U maximum value for phase diagrams
dt    = 0.05;           % Frame duration in seconds (for GIF plotting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init

% Frequency Values
w = linspace(-wcut,wcut,wres); 

% Initial guess for the local Green's function
if MottBIAS
   gloc_0 = 0; % no bath -> no Kondo resonance -> strong Mott bias :)
else
   gloc_0 = phys.BetheLattice(w + 10^(-3)*1i,D); % D is the DOS "radius"
end

%% Single (U,T) point
fprintf('Single point evaluation @ U = %f, T = %f\n\n',U,1/beta)
[gloc,sloc] = DMFT_loop(gloc_0,w,D,U,beta,mloop,mix,err);
Z = phys.Zweight(w,sloc)
I = phys.LuttingerIntegral(w,sloc,gloc)
S = phys.CorrelationStrenght(w,sloc)
if(DoPLOT && DoSPECTRAL)
    [DOS,SELF_ENERGY] = plot.spectral_frame(w,gloc,sloc,U,beta)
end

if Uline
    %% U-driven MIT line [given T]
    fprintf('U-driven span @ T = %f\n\n',1/beta)
    clear('gloc','sloc','Z','I','S')
    i = 0; U = Umin; 
    while U <= Umax 
        i = i + 1;
        fprintf('< U = %f\n',U);
        [gloc{i},sloc{i}] = DMFT_loop(gloc_0,w,D,U,beta,mloop,mix,err,'quiet');
        Z(i) = phys.Zweight(w,sloc{i});
        I(i) = phys.LuttingerIntegral(w,sloc{i},gloc{i});
        S(i) = phys.CorrelationStrenght(w,sloc{i});
        U = U + Ustep;
    end
    if(DoPLOT)
        u_span = plot.Uline(Z,I,beta,Umin,Ustep,Umax)
    end
    if(DoGIF && DoSPECTRAL)
        plot.spectral_gif(w,gloc,sloc,Umin:Ustep:Umax,1/beta,dt)
    end

end

if Tline
    %% T-driven MIT line [given U]
    fprintf('T-driven span @ U = %f\n\n',U)
    clear('gloc','sloc','Z','I','S')
    i = 0; T = Tmin;
    while T <= Tmax 
        i = i + 1; beta = 1/T;
        fprintf('< T = %f\n',T);
        [gloc{i},sloc{i}] = DMFT_loop(gloc_0,w,D,U,beta,mloop,mix,err,'quiet');
        Z(i) = phys.Zweight(w,sloc{i});
        I(i) = phys.LuttingerIntegral(w,sloc{i},gloc{i});
        S(i) = phys.CorrelationStrenght(w,sloc{i});
        T = T + Tstep;
    end
    if(DoPLOT)
        t_span = plot.Tline(Z,U,Tmin,Tstep,Tmax)
    end
    if(DoGIF && DoSPECTRAL)
        plot.spectral_gif(w,gloc,sloc,U,Tmin:Tstep:Tmax,dt)
    end
end

if UTscan
    %% Full Phase-Diagram [U-driven]
    fprintf('Full phase diagram\n\n')
    clear('gloc','sloc','Z','I','S')
    i = 0; T = Tmin; %restart_gloc = gloc_0;
    while T < Tmax
        i = i + 1;
        j = 0; 
        U = Umin; 
        while U <= Umax  
            j = j + 1; beta = 1/T;
            fprintf('< U = %f, T = %f\n',U, T);
            [gloc{i,j},sloc{i,j}] = DMFT_loop(gloc_0,w,D,U,beta,mloop,mix,err,'quiet');
            %restart_gloc = gloc{i,j};
            Z(i,j) = phys.Zweight(w,sloc{i,j});
            I(i,j) = phys.LuttingerIntegral(w,sloc{i,j},gloc{i,j});
            S(i,j) = phys.CorrelationStrenght(w,sloc{i,j});
            U = U + Ustep;
        end
        T = T + Tstep; 
    end
    if(DoPLOT)
        phasemap = plot.phase_diagram(S,Umin,Ustep,Umax,Tmin,Tstep,Tmax)
    end
end
