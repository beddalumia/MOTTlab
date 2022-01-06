function trimFigure(inFig)
%Trims the white borders of the figure when exporting to PDF
if nargin < 1
    inFig = gcf;
end
set(inFig,'PaperPositionMode','auto');
figpos = get(inFig,'PaperPosition');
set(inFig,'PaperSize',[figpos(3) figpos(4)]);
end

