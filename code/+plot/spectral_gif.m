function spectral_gif(w,gloc,sloc,Uvec,Tvec,dt)
%% SPECTRAL_GIF 
%  Builds a GIF for the line-evolution of the spectral functions
%  the function is overloaded to the two types of lines (U- and T-driven)
%  w    : array of frequency values
%  gloc : cell of gloc(w) arrays
%  sloc : cell of sloc(w) arrays
%  Uvec : array of Hubbard-U span values
%  Tvec : array of Temperature span values
%  dt   : delay-time for the GIF frames
%  ------------------------------------------------------------------------
                                                             global DoDEBUG
    fprintf('Start GIF building...\n\n');
    if(length(Uvec)>length(Tvec))
        %% U-driven MIT
        for i = 1:length(Uvec)
            U = Uvec(i);
            beta = 1/Tvec;
            % Build the plots
            [DOS,SE] = plot.spectral_frame(w,gloc{i},sloc{i},U,beta,'invisible');
            % Set filenames
            TitleString = sprintf('beta%d', beta);
            DosName = append('uDOS_',TitleString,'.gif');
            SigmaName = append('uSigma_',TitleString,'.gif');
            % Generate and write the GIF frames 
            plot.push_frame(DosName,i,length(Uvec),dt,DOS);
            plot.push_frame(SigmaName,i,length(Uvec),dt,SE);
            % Close the figures
            close(DOS);
            close(SE);
if DoDEBUG
            [~,LI] = phys.LuttingerIntegral(w,gloc{i},sloc{i});
            title(sprintf('IPT  |  DOS @ U/t = %.2f, beta = %d',2*U,beta));
            LName = append('Luttinger_',TitleString,'.gif');
            plot.push_frame(LName,i,length(Uvec),dt,LI);
            close(LI);
            [~,ZF] = phys.Zweight(w,sloc{i});
            title(sprintf('IPT  |  DOS @ U/t = %.2f, beta = %d',2*U,beta));
            ZName = append('Zfit_',TitleString,'.gif');
            plot.push_frame(ZName,i,length(Uvec),dt,ZF);
            close(ZF);
end
        end
    else
        %% T-driven MIT
        for i = 1:length(Tvec)
            U = Uvec;
            beta = 1/Tvec(i);
            % Build the plots
            [DOS,SE] = plot.spectral_frame(w,gloc{i},sloc{i},U,beta,'invisible');
            % Set filenames
            TitleString = sprintf('U%f', U);
            DosName = append('betaDOS_',TitleString,'.gif');
            SigmaName = append('betaSigma_',TitleString,'.gif');
            % Generate and write the GIF frames 
            plot.push_frame(DosName,i,length(Tvec),dt,DOS);
            plot.push_frame(SigmaName,i,length(Tvec),dt,SE);
            % Close the figures
            close(DOS);
            close(SE);
        end
    end
    fprintf('...GIFs have been built.\n\n');
end 
   
    
