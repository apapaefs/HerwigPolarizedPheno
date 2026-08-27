BEGIN PLOT /COMPASS_2017_I1444985/.*
LogY=0
RatioPlot=0
XLabel=$z$
YLabel=$dM^h(x,y,z)/dz$
LegendAlign=r
END PLOT

BEGIN PLOT /COMPASS_2017_I1444985/Multiplicity_.*_cells
XLabel=Published sparse $(x,y,z)$ cell
END PLOT

BEGIN PLOT /COMPASS_2017_I1444985/Pull_.*
XLabel=Published sparse cell
YLabel=$(M^{\mathrm{theory}}-M^{\mathrm{data}})/\sigma$
YMin=-5
YMax=5
ConnectBins=0
END PLOT

BEGIN PLOT /COMPASS_2017_I1444985/DIAGNOSTICS/.*
XLabel=Diagnostic bin
YLabel=Generator diagnostic
END PLOT
