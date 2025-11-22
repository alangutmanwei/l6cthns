function figQuality(figh, axh, dim)
    %figh = figure handle
    %axh = axis handle
    %dim = [w h] in inches
    
    figh.Units = 'inches';
    figh.Position = [0 0 dim];
    figh.Renderer = 'Painters';
    figh.Color = [1 1 1];
    figh.PaperPositionMode = 'auto'; 
    box off;
    %set(axh,'LooseInset',get(axh,'TightInset')) 
    
    axh.XColor = 'k'; axh.YColor = 'k'; axh.ZColor = 'k';
    allAxesInFigure = findall(figh,'type','axes');
    for i=1:size(allAxesInFigure,1)
        axh = allAxesInFigure(i);
        axh.FontSizeMode = 'manual';
        axh.FontSize = 12;
        plblsz = 16; %preferred label size
        axh.LabelFontSizeMultiplier = plblsz / axh.FontSize;
        axh.FontName = 'Helvetica';
        axh.LineWidth = 1.5;
        axh.TickDir = 'out';
    end
end
