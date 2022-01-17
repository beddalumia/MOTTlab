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
            % Build the plots
            [DOS,SE] = plot.spectral_frame(w,gloc{i},sloc{i},U,beta,'invisible');
            % Set filenames
            TitleString = sprintf('beta%d', beta);
            DosName = append('uDOS_',TitleString,'.gif');
            SigmaName = append('uSigma_',TitleString,'.gif');
            % Generate and write the GIF frames 
            push_gif_frame(DosName,i,length(Uvec),dt,DOS);
            push_gif_frame(SigmaName,i,length(Uvec),dt,SE);
            % Close the figures
            close(DOS);
            close(SE);
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
            push_gif_frame(DosName,i,length(Tvec),dt,DOS);
            push_gif_frame(SigmaName,i,length(Tvec),dt,SE);
            % Close the figures
            close(DOS);
            close(SE);
        end
    end
end 
   
function push_gif_frame(gifname,iframe,nframe,dt,infig)
% gifname   :   filename for the gif, should end with .gif
% iframe    :   id number of the given frame, as an integer
% nframe    :   total number of expected frames, as an integer
% dt        :   duration for the given frame, in seconds
% infig     :   optional handle for the figure to print
                                                        if nargin < 4
                                                           infig = gcf;
                                                        end        
     % Capture the plot as an image
     im = print(infig,'-RGBImage');
     % Suitable conversion
     [ind,cm] = rgb2ind(im,256,'nodither');
     % Write to the GIF file
     if iframe == 1
        imwrite(ind,cm,gifname,'gif','Loopcount',inf,'DelayTime',dt);
     else
        imwrite(ind,cm,gifname,'gif','WriteMode','append','DelayTime',dt);
     end
     % Write info to stdout
     info = sprintf('Added %d-th frame of %d to ',iframe,nframe);
     info = [info,'<',gifname,'>','\n'];
     fprintf(info);
end
   
    
