function [DOS,SELF_ENERGY] = spectral_frame(w,gloc,sloc,U,beta,fmode)
%% SPECTRAL_FRAME
%  Builds nice plots for the spectral functions
%  w     : array of frequency values [float]
%  gloc  : array of G_loc spectral values [complex]
%  sloc  : array of Sigma_loc spectral values [complex]
%  U     : single Hubbard-U value [float]
%  beta  : single inverse-temperature value [float]
%  fmode : string flag to set figures visibility [OPTIONAL, default: 'show']
%  ------------------------------------------------------------------------
    if(~exist('fmode','var'))
        fmode = 'show';
    end
    show = strcmp(fmode,'show');
    DOS = figure("Name",'Spectral Function (DOS)','Visible','off');
    if(show)
        set(DOS,'Visible', 'on');
    end
    FilledStates = w(w<=0);
    FilledDOS = -imag(gloc(w<=0))/pi;
    EmptyStates = w(w>0);
    EmptyDOS = -imag(gloc(w>0))/pi; 
    area(FilledStates, FilledDOS, 'FaceColor', [0.7 0.7 0.7]); hold on
    area(EmptyStates, EmptyDOS, 'FaceColor', [1 1 1]);
    title(sprintf('IPT  |  DOS @ U/D = %.2f, beta = %d', U, beta));
    xlabel('\omega')
    ylabel('A(\omega)')
    ylim([0,1])
    legend('Filled States','Empty States')
    SELF_ENERGY = figure("Name",'Self-Energy','Visible','off');
    if(show)
        set(SELF_ENERGY,'Visible', 'on');
    end
    plot(w,real(sloc)./U,'LineWidth',1); hold on
    plot(w,imag(sloc)./U,'LineWidth',1);
    title(sprintf('IPT  |  Self-Energy @ U/D = %.2f, beta = %d', U, beta)); 
    xlabel('\omega')
    ylabel('Units of U')
    ylim([-1,1])
    legend('Re \Sigma_{loc}(\omega)','Im \Sigma_{loc}(\omega)')
end
