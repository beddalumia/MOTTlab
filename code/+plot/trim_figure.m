function trim_figure(infig)
%Trims the white borders of the figure when exporting to PDF
if nargin < 1
    infig = gcf;
end
set(infig,'PaperPositionMode','auto');
figpos = get(infig,'PaperPosition');
set(infig,'PaperSize',[figpos(3) figpos(4)]);
end

