BEGIN PLOT /COMPASS_2017_I1483098/.*
LogY=0
RatioPlot=0
XLabel=$z$
YLabel=$dM^K(x,y,z)/dz$
LegendAlign=r
END PLOT

BEGIN PLOT /COMPASS_2017_I1483098/Multiplicity_.*_cells
XLabel=Published sparse $(x,y,z)$ cell
END PLOT

BEGIN PLOT /COMPASS_2017_I1483098/Pull_.*
XLabel=Published sparse cell
YLabel=$(M^{\mathrm{theory}}-M^{\mathrm{data}})/\sigma$
YMin=-5
YMax=5
ConnectBins=0
END PLOT

BEGIN PLOT /COMPASS_2017_I1483098/DIAGNOSTICS/.*
XLabel=Diagnostic bin
YLabel=Generator diagnostic
END PLOT
