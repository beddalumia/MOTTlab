function spectral_gif(w,gloc,sloc,Uvec,Tvec,dt)
%% Builds a GIF for the line-evolution of the spectral functions
%  the function is overloaded to the two types of lines (U- and T-driven)
%  w    : array of frequency values
%  gloc : cell of gloc(w) arrays
%  sloc : cell of sloc(w) arrays
%  Uvec : array of Hubbard-U span values
%  Tvec : array of Temperature span values
%  dt   : delay-time for the GIF frames
%  ------------------------------------------------------------------------
    fprintf('Start GIF building...\n\n');
    if(length(Uvec)>length(Tvec))
        %% U-driven MIT
        for i = 1:length(Uvec)
            U = Uvec(i);
            beta = 1/Tvec;
            % Build the plot
            [DOS,SE] = plot.spectral_frame(w,gloc{i},sloc{i},U,beta,false);
            % Capture the plot as an image 
            imDOS = print(DOS,'-RGBImage');
            imSE = print(SE,'-RGBImage'); 
            [imDOSind,cmDOS] = rgb2ind(imDOS,256,'nodither');
            [imSEind,cmSE] = rgb2ind(imSE,256,'nodither');
            % Set filenames
            TitleString = sprintf('beta%d', beta);
            DosName = append('uDOS_',TitleString,'.gif');
            SigmaName = append('uSigma_',TitleString,'.gif');
            % Write to the GIF files 
            if i == 1 
               imwrite(imDOSind,cmDOS,DosName,'gif', 'Loopcount',inf,'DelayTime',dt);
               imwrite(imSEind,cmSE,SigmaName,'gif', 'Loopcount',inf,'DelayTime',dt);
            else 
               imwrite(imDOSind,cmDOS,DosName,'gif','WriteMode','append','DelayTime',dt);
               imwrite(imSEind,cmSE,SigmaName,'gif','WriteMode','append','DelayTime',dt);
            end
            fprintf('Added %d-th frame of %d\n',i,length(Uvec));
            close(DOS);
            close(SE);
        end
    else
        %% T-driven MIT
        for i = 1:length(Tvec)
            U = Uvec;
            beta = 1/Tvec(i);
            % Build the plot
            [DOS,SE] = plot.spectral_frame(w,gloc{i},sloc{i},U,beta,false);
            % Capture the plot as an image 
            imDOS = print(DOS,'-RGBImage');
            imSE = print(SE,'-RGBImage'); 
            [imDOSind,cmDOS] = rgb2ind(imDOS,256,'nodither');
            [imSEind,cmSE] = rgb2ind(imSE,256,'nodither');
            % Set filenames
            TitleString = sprintf('U%f', U);
            DosName = append('betaDOS_',TitleString,'.gif');
            SigmaName = append('betaSigma_',TitleString,'.gif');
            % Write to the GIF files 
            if i == 1 
               imwrite(imDOSind,cmDOS,DosName,'gif', 'Loopcount',inf,'DelayTime',dt);
               imwrite(imSEind,cmSE,SigmaName,'gif', 'Loopcount',inf,'DelayTime',dt);
            else 
               imwrite(imDOSind,cmDOS,DosName,'gif','WriteMode','append','DelayTime',dt);
               imwrite(imSEind,cmSE,SigmaName,'gif','WriteMode','append','DelayTime',dt);
            end 
            fprintf('Added %d-th frame of %d\n',i,length(Tvec));
            close(DOS);
            close(SE);
        end
    end
end 
   
   
   
    
