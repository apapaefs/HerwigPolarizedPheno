BEGIN PLOT /STAR_2021_I1850855/*
LogY=0
RatioPlot=0
YLabel=$A_{LL}$
LegendAlign=l
LegendAnchor=upper left
LegendXPos=0.00
LegendYPos=0.95
END PLOT

BEGIN PLOT /STAR_2021_I1850855/Pull_*
Title=STAR data--prediction pull
YLabel=$(A_{LL}^{\mathrm{theory}}-A_{LL}^{\mathrm{data}})/\sigma$
ConnectBins=0
YMin=-5
YMax=5
END PLOT

BEGIN PLOT /STAR_2021_I1850855/inclusive_*
XLabel=Parton-jet $p_T$ [GeV]
END PLOT

BEGIN PLOT /STAR_2021_I1850855/dijet_*
XLabel=Parton-dijet mass [GeV]
END PLOT
