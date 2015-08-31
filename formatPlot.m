function [] = formatPlot(varargin)
% Run formatPlot after all plots have been created in order to format every
% open, visible plot using the given style. If no input argument is given,
% the default style will be used.
%
% INPUTS:
% style = 'paper'
%
% Rick Fenrich 8/31/15

if(nargin == 0) % default
    fontName = 'CMU Serif';
    titleFontSize = 16;
    axisLabelFontSize = 14;
    axisFontSize = 12;
    legendFontSize = 12;
elseif(strcmp(varargin,'paper'))
    fontName = 'CMU Serif';
    titleFontSize = 16;
    axisLabelFontSize = 14;
    axisFontSize = 12;
    legendFontSize = 12;
else
    error('Style not defined.');
end

allFigures = get(0,'children'); % obtains array of all figures from root

for ii = 1:length(allFigures) % for each figure fix the labels
    figure(ii)
    fixLabels();
end

function [] = fixLabels()
    
    currentFigChildren = get(gcf,'children');
    
    for jj = 1:length(currentFigChildren)
        
        % Format axis
        if(strcmp(currentFigChildren(jj).Type,'axes')) % if object is an axis

            % Format title
            set(currentFigChildren(jj).Title,'FontName',fontName);
            set(currentFigChildren(jj).Title,'FontSize',titleFontSize);
            set(currentFigChildren(jj).Title,'FontWeight','bold');

            % Format axes
            set(currentFigChildren(jj),'FontName',fontName);
            set(currentFigChildren(jj),'FontSize',axisFontSize);
            set(currentFigChildren(jj).XLabel,'FontSize',axisLabelFontSize);
            set(currentFigChildren(jj).YLabel,'FontSize',axisLabelFontSize);
            
        % Format legend
        elseif(strcmp(currentFigChildren(jj).Type,'legend')) % if object is a legend

            % Format legend
            set(currentFigChildren(jj),'FontSize',legendFontSize);
            
        end
        
    end

end

end

