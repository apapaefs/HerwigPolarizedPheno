BEGIN PLOT /COMPASS_2009_I820721/A1_.*
LogX=1
LogY=0
RatioPlot=0
XLabel=$x$
YLabel=$A_{1,d}^{h}$
LegendAlign=l
END PLOT

BEGIN PLOT /COMPASS_2009_I820721/Pull_.*
LogX=1
RatioPlot=0
XLabel=$x$
YLabel=$(A_1^{\mathrm{theory}}-A_1^{\mathrm{data}})/\sigma$
YMin=-5
YMax=5
ConnectBins=0
END PLOT

BEGIN PLOT /COMPASS_2009_I820721/DIAGNOSTICS/.*
RatioPlot=0
XLabel=Diagnostic bin
YLabel=Generator diagnostic
END PLOT
