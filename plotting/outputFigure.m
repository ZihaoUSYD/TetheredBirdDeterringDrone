function outputFigure(figureFileName,figureWidthScale,figureHeightScale)

latexFontSize = 12;
latexTextwidth = 6.53278;
latex_fig(latexFontSize,figureWidthScale*latexTextwidth,figureHeightScale*latexTextwidth);
export_fig(figureFileName,'-pdf');

end