%% BSD 3-Clause License
% 
% Copyright (c) 2020-2022, Gabriele Bellomia
% All rights reserved.

clearvars; clc;
global DEBUG FAST

try
    pkg load signal     % GNU Octave option
end

%% INPUT: Physical Parameters 
U    = 4.0;             % On-site Repulsion
beta = 1e3;             % Inverse Temperature
D    = 1.0;             % Noninteracting half-bandwidth
latt = 'bethe';         % Noninteracting band-dispersion 
                        % ['bethe','cubic','square','chain'...]

%% INPUT: Boolean Flags
MottBIAS     = 0;       % Changes initial guess of gloc (strongly favours Mott phase)
ULINE        = 1;       % Takes and fixes the given beta value and performs a U-driven line
TLINE        = 0;       % Takes and fixes the given U value and performs a T-driven line
UTSCAN       = 0;       % Ignores both given U and beta values and builds a full phase diagram
SPECTRAL     = 0;       % Controls plotting of spectral functions
PLOT         = 1;       % Controls plotting of *all static* figures
GIF          = 0;       % Controls plotting of *animated* figures
PRINT        = 0;       % Controls file printing (for single points)
UARRAY       = 0;       % Activates SLURM scaling of interaction values
TARRAY       = 0;       % Activates SLURM scaling of temperature values 
RESTART      = 1;       % Activates the restarting strategies for lines               
DEBUG        = 0;       % Activates debug prints / plots / operations
FAST         = 1;       % Activates fast FFTW-based convolutions

%% INPUT: Control Parameters
mloop = 1000;           % Max number of DMFT iterations 
err   = 1e-5;           % Convergence threshold for self-consistency
mix   = 0.30;           % Mixing parameter for DMFT iterations (=1 means full update)
wres  = 2^15;           % Energy resolution in real-frequency axis
wcut  = 6.00;           % Energy cutoff in real-frequency axis
vcut  = 6.00;           % Energy cutoff in imag-frequency axis
Umin  = 0.00;           % Hubbard U minimum value for phase diagrams
Ustep = 0.10;           % Hubbard U incremental step for phase diagrams
Umax  = 4.50;           % Hubbard U maximum value for phase diagrams
Tmin  = 1e-3;           % Temperature U minimum value for phase diagrams
Tstep = 1e-3;           % Temperature incremental step for phase diagrams
Tmax  = 5e-2;           % Temperature U maximum value for phase diagrams
dt    = 0.05;           % Frame duration in seconds (for GIF plotting)

%% SLURM-SCALING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aID = str2double(getenv('SLURM_ARRAY_TASK_ID'));

if UARRAY
   U = U*aID;
end

if TARRAY
   beta = beta*aID;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init

% Frequency Values
w = linspace(-wcut,wcut,wres); 

% Initial guess for the local Green's function
if MottBIAS
   seed = 0; % no bath -> no Kondo resonance -> strong Mott bias :)
else
   seed = phys.gloc(w + 10^(-3)*1i,D,latt);
end

%% Workflows

if not( ULINE || TLINE || UTSCAN )
    %% Single (U,T) point
    fprintf('Single point evaluation @ U = %f, T = %f\n\n',U,1/beta); tic
    [gloc,sloc] = dmft_loop(seed,w,U,beta,D,latt,mloop,mix,err);
    [iv,gmatsu] = phys.matsubara(w,gloc,beta,vcut); m = imag(gmatsu(1));
    [iv,smatsu] = phys.matsubara(w,sloc,beta,vcut);
    Z = phys.zetaweight(w,sloc); 
    z = phys.zetaweight(iv,smatsu);
    I = phys.luttinger(w,sloc,gloc);
    S = phys.strcorrel(w,sloc);
    if(PLOT && SPECTRAL)
        [DOS,SELF_ENERGY] = plot.spectral_frame(w,gloc,sloc,U,beta,D);
    end
    ET = [0,0,toc]; fmt = 'hh:mm:ss.SSS';
    fprintf('> %s < elapsed time\n\n',duration(ET,'format',fmt));
    if PRINT
        writematrix(U,sprintf('U%f_HU.dat',U));
        writematrix(I,sprintf('U%f_IL.dat',U)); 
        writematrix(Z,sprintf('U%f_ZF.dat',U)); 
        writematrix(S,sprintf('U%f_SR.dat',U));
        writematrix(m,sprintf('U%f_MT.dat',U));
    end
end

if ULINE
    %% U-driven MIT line [given T]
    fprintf('U-driven span @ T = %f\n\n',1/beta); tic
    clear('gloc','sloc','Z','I','S'); 
    Uvec = Umin:Ustep:Umax; NU = length(Uvec);
    gloc_0 = seed; gloc = cell(NU,1); sloc = gloc; gmatsu = gloc; smatsu = sloc;
    Z = zeros(NU,1); I = zeros(NU,1); S = zeros(NU,1); m = Z; z = Z; d = Z; n = Z; d = Z;
    for i = 1:NU 
        U = Uvec(i);
        fprintf('< U = %f\n',U);
        [gloc{i},sloc{i}] = dmft_loop(gloc_0,w,U,beta,D,latt,mloop,mix,err,'quiet');
        [iv,gmatsu{i}] = phys.matsubara(w,gloc{i},beta,vcut); m(i) = imag(gmatsu{i}(1));
        [iv,smatsu{i}] = phys.matsubara(w,sloc{i},beta,vcut);
        if(RESTART)
           gloc_0 = gloc{i}; 
        end
        Z(i) = phys.zetaweight(w,sloc{i}); 
        z(i) = phys.zetaweight(iv,smatsu{i});
        d(i) = phys.docc(iv,smatsu{i},gmatsu{i},U);
        I(i) = phys.luttinger(w,sloc{i},gloc{i});
        S(i) = phys.strcorrel(w,sloc{i});
    end
    if(PLOT)
        u_span = plot.Uline(Z,d,beta,Umin,Ustep,Umax,D);
    end
    if(GIF && SPECTRAL)
        plot.spectral_gif(w,gloc,sloc,Umin:Ustep:Umax,1/beta,D,dt);
    end
    ET = [0,0,toc]; fmt = 'hh:mm:ss.SSS';
    fprintf('> %s < elapsed time\n\n',duration(ET,'format',fmt));
end

if TLINE
    %% T-driven MIT line [given U]
    fprintf('T-driven span @ U = %f\n\n',U); tic
    clear('gloc','sloc','Z','I','S')
    Tvec = Tmin:Tstep:Tmax; NT = length(Tvec);
    gloc_0 = seed; gloc = cell(NT,1); sloc = gloc; gmatsu = gloc;
    Z = zeros(NT,1); I = zeros(NT,1); S = zeros(NT,1); m = zeros(NT,1);
    for i = 1:NT 
        T = Tvec(i); beta = 1/T;
        fprintf('< T = %f\n',T);
        [gloc{i},sloc{i}] = dmft_loop(gloc_0,w,U,beta,D,latt,mloop,mix,err,'quiet');
        [iv,gmatsu{i}] = phys.matsubara(w,gloc{i},beta,vcut); m(i) = imag(gmatsu{i}(1));
        if(RESTART)
           gloc_0 = gloc{i}; 
        end
        Z(i) = phys.zetaweight(w,sloc{i});
        I(i) = phys.luttinger(w,sloc{i},gloc{i});
        S(i) = phys.strcorrel(w,sloc{i});
    end
    if(PLOT)
        t_span = plot.Tline(Z,U,Tmin,Tstep,Tmax,D);
    end
    if(GIF && SPECTRAL)
        plot.spectral_gif(w,gloc,sloc,U,Tmin:Tstep:Tmax,D,dt);
    end
    ET = [0,0,toc]; fmt = 'hh:mm:ss.SSS';
    fprintf('> %s < elapsed time\n\n',duration(ET,'format',fmt));
end

if UTSCAN
    %% Full Phase-Diagram [U-driven]
    fprintf('Full phase diagram\n\n'); tic; 
    clear('gloc','sloc','Z','I','S')
    Tvec = Tmin:Tstep:Tmax; NT = length(Tvec);
    Uvec = Umin:Ustep:Umax; NU = length(Uvec);
    gloc = cell(NT,NU); sloc = gloc; gmatsu = gloc;
    Z = zeros(NT,NU);  I = Z; S = Z; m = Z;
    feature('numcores');
    parfor i = 1:NT 
        T = Tvec(i); beta = 1/T;
        gloc_0 = seed;
        for j = 1:NU  
            U = Uvec(j);
            fprintf('< U = %f, T = %f\n',U, T);
            [gloc{i,j},sloc{i,j}] = dmft_loop(gloc_0,w,U,beta,D,latt,mloop,mix,err,'quiet');
            [iv,gmatsu{i,j}] = phys.matsubara(w,gloc{i,j},beta,vcut); m(i,j) = imag(gmatsu{i,j}(1));
            if(RESTART)
               gloc_0 = gloc{i,j};
            end
            Z(i,j) = phys.zetaweight(w,sloc{i,j});
            I(i,j) = phys.luttinger(w,sloc{i,j},gloc{i,j});
            S(i,j) = phys.strcorrel(w,sloc{i,j});
        end
    end
    if(PLOT)
        S = S/max(max(S));
        zetamap = plot.phase_diagram(Z,Umin,Ustep,Umax,Tmin,Tstep,Tmax,D);
        luttmap = plot.phase_diagram(I,Umin,Ustep,Umax,Tmin,Tstep,Tmax,D);
        strcmap = plot.phase_diagram(S,Umin,Ustep,Umax,Tmin,Tstep,Tmax,D);
    end
    ET = [0,0,toc]; fmt = 'hh:mm:ss.SSS';
    fprintf('> %s < elapsed time\n\n',duration(ET,'format',fmt));
end
