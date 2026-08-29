BEGIN PLOT /MC_POLDIJETS/.*
RatioPlot=0
LogY=0
ConnectBins=1
LegendAlign=l
LegendAnchor=upper left
LegendXPos=0.00
LegendYPos=0.95
END PLOT

BEGIN PLOT /MC_POLDIJETS/ALL_.*
YLabel=$A_{LL}$
END PLOT

BEGIN PLOT /MC_POLDIJETS/DeltaSigmaLL_.*
YLabel=$\mathrm{d}\Delta\sigma_{LL}/\mathrm{d}X$ [pb / unit]
END PLOT

BEGIN PLOT /MC_POLDIJETS/SigmaUU_.*
YLabel=$\mathrm{d}\sigma_{UU}/\mathrm{d}X$ [pb / unit]
LogY=1
END PLOT

BEGIN PLOT /MC_POLDIJETS/SigmaPP_.*
YLabel=$\mathrm{d}\sigma^{++}/\mathrm{d}X$ [pb / unit]
LogY=1
END PLOT

BEGIN PLOT /MC_POLDIJETS/SigmaPM_.*
YLabel=$\mathrm{d}\sigma^{+-}/\mathrm{d}X$ [pb / unit]
LogY=1
END PLOT

BEGIN PLOT /MC_POLDIJETS/SigmaMP_.*
YLabel=$\mathrm{d}\sigma^{-+}/\mathrm{d}X$ [pb / unit]
LogY=1
END PLOT

BEGIN PLOT /MC_POLDIJETS/SigmaMM_.*
YLabel=$\mathrm{d}\sigma^{--}/\mathrm{d}X$ [pb / unit]
LogY=1
END PLOT

BEGIN PLOT /MC_POLDIJETS/ALL_dijet_rate$
XLabel=Selected loose-dijet sample
YLabel=$A_{LL}$
LogY=0
END PLOT

BEGIN PLOT /MC_POLDIJETS/DeltaSigmaLL_dijet_rate$
XLabel=Selected loose-dijet sample
YLabel=$\Delta\sigma_{LL}$ [pb]
LogY=0
END PLOT

BEGIN PLOT /MC_POLDIJETS/SigmaUU_dijet_rate$
XLabel=Selected loose-dijet sample
YLabel=$\sigma_{UU}$ [pb]
LogY=0
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet1_pt$
XLabel=Leading parton-jet $p_T$ [GeV]
XMax=70
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet2_pt$
XLabel=Subleading parton-jet $p_T$ [GeV]
XMax=45
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet3_pt$
XLabel=Third parton-jet $p_T$ [GeV]
XMax=20
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet4_pt$
XLabel=Fourth parton-jet $p_T$ [GeV]
XMax=15
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet1_eta$
XLabel=Leading parton-jet $\eta$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet2_eta$
XLabel=Subleading parton-jet $\eta$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet3_eta$
XLabel=Third parton-jet $\eta$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet4_eta$
XLabel=Fourth parton-jet $\eta$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_pt_average$
XLabel=$(p_{T,1}+p_{T,2})/2$ [GeV]
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_pt_ratio_21$
XLabel=$p_{T,2}/p_{T,1}$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_pt_ratio_31$
XLabel=$p_{T,3}/p_{T,1}$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_dijet_mass$
XLabel=Parton-dijet mass [GeV]
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_dijet_pt$
XLabel=Parton-dijet $p_T$ [GeV]
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_dijet_ht$
XLabel=$p_{T,1}+p_{T,2}$ [GeV]
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_delta_phi$
XLabel=$\Delta\phi(j_1,j_2)$ [rad]
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_abs_delta_eta$
XLabel=$|\Delta\eta(j_1,j_2)|$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_delta_r$
XLabel=$\Delta R(j_1,j_2)$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_delta_r_13$
XLabel=$\Delta R(j_1,j_3)$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_delta_r_24$
XLabel=$\Delta R(j_2,j_4)$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_delta_r_14$
XLabel=$\Delta R(j_1,j_4)$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_eta_boost$
XLabel=$(\eta_1+\eta_2)/2$
END PLOT

BEGIN PLOT /MC_POLDIJETS/(ALL|DeltaSigmaLL|SigmaUU)_cos2_delta_phi$
XLabel=$\cos(2\Delta\phi(j_1,j_2))$
END PLOT

BEGIN PLOT /MC_POLDIJETS/DIAGNOSTICS/SingleSpinA_.*
YLabel=$A_L^{(1)}$
END PLOT

BEGIN PLOT /MC_POLDIJETS/DIAGNOSTICS/SingleSpinB_.*
YLabel=$A_L^{(2)}$
END PLOT

BEGIN PLOT /MC_POLDIJETS/DIAGNOSTICS/Parity_.*
YLabel=Parity residual
END PLOT
