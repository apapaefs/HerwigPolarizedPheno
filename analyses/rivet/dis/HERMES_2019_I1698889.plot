BEGIN PLOT /HERMES_2019_I1698889/.*
LogY=0
RatioPlot=0
YLabel=$A_{\parallel}^{h}$
LegendAlign=l
LegendAnchor=upper left
LegendXPos=0.05
LegendYPos=0.95
END PLOT

BEGIN PLOT /HERMES_2019_I1698889/Pull_.*
Title=HERMES data--prediction pull
YLabel=$(A_{\parallel}^{\mathrm{theory}}-A_{\parallel}^{\mathrm{data}})/\sigma$
ConnectBins=0
YMin=-5
YMax=5
END PLOT

BEGIN PLOT /HERMES_2019_I1698889/.*_x$
XLabel=$x$
END PLOT

BEGIN PLOT /HERMES_2019_I1698889/.*_xz$
XLabel=Flattened $(x,z)$ bin
END PLOT

BEGIN PLOT /HERMES_2019_I1698889/.*_xpt$
XLabel=Flattened $(x,P_{hT})$ bin
END PLOT

BEGIN PLOT /HERMES_2019_I1698889/.*_xzpt$
XLabel=Flattened $(x,z,P_{hT})$ bin
END PLOT

BEGIN PLOT /HERMES_2019_I1698889/DIAGNOSTICS/.*
YLabel=Generator diagnostic
XLabel=Diagnostic bin
END PLOT

BEGIN PLOT /HERMES_2019_I1698889/DIAGNOSTICS/UnpolarizedCosPhi_.*
YLabel=$2\langle\cos\phi\rangle_{UU}$
XLabel=Unpolarized yield-moment bin
END PLOT

BEGIN PLOT /HERMES_2019_I1698889/PUBLISHED_AParallelCosPhi/.*
YLabel=$A_{\parallel}^{h,\cos\phi}$
XLabel=Published spin-asymmetry fit bin
END PLOT
