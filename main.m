%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

clear all; clc

%% INPUT: Physical Parameters 
D    = 1;               % Bandwidth
U    = 3;               % On-site Repulsion    } Overriden if PhaseDiagram
beta = 10^3;            % Inverse Temperature  } flag is set to true...

%% INPUT: Boolean Flags
MottBIAS     = 0;       % Changes initial guess of gloc (strongly favours Mott phase)
Uline        = 0;       % Takes and fixes the given beta value and performs a U-driven line
Tline        = 0;       % Takes and fixes the given U value and performs a T-driven line
UTscan       = 1;       % Ignores both given U and beta values and builds a full phase diagram
DoSPECTRAL   = 0;       % Controls plotting of spectral functions
DoPLOT       = 1;       % Controls plotting of *all* the figures

%% INPUT: Control Parameters
mloop = 1000;           % Max number of DMFT iterations 
err   = 10^(-1);        % Convergence threshold for self-consistency
mix   = 0.1;            % Mixing parameter for DMFT iterations (=1 means full update)
wres  = 2^12;           % Energy resolution in real-frequency axis
Ustep = 0.05;           % Hubbard U incremental step for phase diagrams
Tstep = 0.005;          % Temperature incremental step for phase diagrams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init

% Frequency Values
w = linspace(-6, 6, wres); 

% Initial guess for the local Green's function
if MottBIAS
   gloc_0 = 0; % no bath -> no Kondo resonance -> strong Mott bias :)
else
   gloc_0 = BetheHilbert(w + 10^(-3)*1i,D); % D is the DOS "radius"
end

%% Single (U,T) point
[gloc,Sigma_loc] = DMFT_loop(gloc_0,w,D,U,beta,mloop,mix,err,false);
Z = Zweight(w,Sigma_loc)
I = LuttingerIntegral(w,Sigma_loc,gloc)
if(DoPLOT && DoSPECTRAL)
    plot_spectral(w,gloc,Sigma_loc,U,beta)
end

if Uline
    %% U-driven MIT line [given T]
    clear('gloc','Sigma_loc','Z','I')
    i = 0; U = 0; 
    while U <= 5.00 
        i = i + 1;
        [gloc{i},Sigma_loc{i}] = DMFT_loop(gloc_0,w,D,U,beta,mloop,mix,err,true);
        Z(i) = Zweight(w,Sigma_loc{i});
        I(i) = LuttingerIntegral(w,Sigma_loc{i},gloc{i});
        U = U + Ustep;
    end
    if(DoPLOT)
        plot_Uline(Z,I,beta,Ustep)
    end
end

if Tline
    %% T-driven MIT line [given U]
    clear('gloc','Sigma_loc','Z','I')
    i = 0; T = 10^(-6);
    while T <= 0.05 
        i = i + 1; beta = 1/T;
        [gloc{i},Sigma_loc{i}] = DMFT_loop(gloc_0,w,D,U,beta,mloop,mix,err,true);
        Z(i) = Zweight(w,Sigma_loc{i});
        I(i) = LuttingerIntegral(w,Sigma_loc{i},gloc{i});
        T = T + Tstep;
    end
    if(DoPLOT)
        plot_Tline(Z,I,U,Tstep)
    end
end

if UTscan
    %% Full Phase-Diagram [U-driven]
    clear('gloc','Sigma_loc','Z','I')
    i = 0; T = 10^(-6); %restart_gloc = gloc_0;
    while T < 0.05
        i = i + 1;
        j = 0; 
        U = 0; 
        while U <= 5.00  
            j = j + 1; beta = 1/T;
            [gloc{i,j},Sigma_loc{i,j}] = DMFT_loop(gloc_0,w,D,U,beta,mloop,mix,err,true);
            %restart_gloc = gloc{i,j};
            Z(i,j) = Zweight(w,Sigma_loc{i,j});
            I(i,j) = LuttingerIntegral(w,Sigma_loc{i,j},gloc{i,j});
            U = U + Ustep;
        end
        T = T + Tstep; 
    end
    if(DoPLOT)
        plot_phase_diagram(Z,Ustep,Tstep)
    end
end

%% Plotter Module
%  -> to be moved in a proper MATLAB package 

function plot_Uline(Z,I,beta,Ustep)
    figure("Name",'U-driven MIT')
    i=0; U = 0;
    while U <= 5.00 
        i = i + 1;
        scatter(2*U,Z(i),'k','filled'); hold on % (Units: D=2t)
        U = U + Ustep;
    end
    title(sprintf('IPT  |  Quasiparticle weight at T = %f',1/beta))
    xlabel('$U/t$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
end

function plot_Tline(Z,I,U,Ustep)
    figure("Name",'T-driven MIT')
    i=0; T = 10^(-6);
    while T <= 0.05
        i = i + 1;
        scatter(T,Z(i),'k','filled'); hold on
        T = T + Ustep;
    end
    title(sprintf('IPT  |  Quasiparticle weight at U/t = %f',2*U)) % (Units: D=2t)
    xlabel('$T$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
end

function plot_phase_diagram(order_parameter,Ustep,Tstep)
    figure("Name",'Phase Diagram')
    X = 0:Ustep:5;
    Y = 10^(-6):Tstep:0.05;
    Z = order_parameter;
    % Remove machine zeros
    Z(Z<10^(-5)) = 0;
    % Surface plot
    surf(X,Y,Z);
end

function plot_spectral(w,gloc,sloc,U,beta)
    figure("Name",'Spectral Function (DOS)')
    FilledStates = w(w<=0);
    FilledDOS = -imag(gloc(w<=0));
    EmptyStates = w(w>0);
    EmptyDOS = -imag(gloc(w>0)); 
    area(FilledStates, FilledDOS, 'FaceColor', [0.7 0.7 0.7]); hold on
    area(EmptyStates, EmptyDOS, 'FaceColor', [1 1 1]);
    title(sprintf('IPT  |  DOS @ U/t = %.2f, beta = %d', 2*U, beta)); % (Units: D=2t)
    xlabel('\omega')
    ylabel('\pi A(\omega)')
    legend('Filled States','Empty States')
    figure("Name",'Self-Energy')
    plot(w,real(sloc),'LineWidth',1); hold on
    plot(w,imag(sloc),'LineWidth',1);
    title(sprintf('IPT  |  Self-Energy @ U/t = %.2f, beta = %d', 2*U, beta)); % (Units: D=2t)
    xlabel('\omega')
    legend('Re \Sigma_{loc}(\omega)','Im \Sigma_{loc}(\omega)')
end
