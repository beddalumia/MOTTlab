function push_frame(gifname,iframe,nframe,dt,infig)
%% PUSH_FRAME: appends a frame to a .gif file [high quality]
% gifname   :   filename for the gif, should end with .gif
% iframe    :   id number of the given frame, as an integer
% nframe    :   total number of expected frames, as an integer
% dt        :   duration for the given frame, in seconds
% infig     :   optional handle for the figure to print
                                                        if nargin < 5
                                                           infig = gcf;
                                                        end        
     % Capture the plot as an image
     im = print(infig,'-RGBImage');
     % Suitable conversion
     [ind,cm] = rgb2ind(im,256,'nodither');
     % Prepare output path
     if ~isfolder('../output'),mkdir('../output');end
     path = ['../output/',gifname];
     % Write to the GIF file
     if iframe == 1
        imwrite(ind,cm,path,'gif','Loopcount',inf,'DelayTime',dt);
     else
        imwrite(ind,cm,path,'gif','WriteMode','append','DelayTime',dt);
     end
     % Write info to stdout
     info = sprintf('Added %d-th frame of %d to ',iframe,nframe);
     info = [info,'<',gifname,'>','\n'];
     fprintf(info);
end