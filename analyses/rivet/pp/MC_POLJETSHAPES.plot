BEGIN PLOT /MC_POLJETSHAPES/.*
RatioPlot=0
LogY=0
ConnectBins=1
LegendAlign=l
LegendAnchor=upper left
LegendXPos=0.00
LegendYPos=0.95
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/ALL_.*
YLabel=$A_{LL}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/DeltaSigmaLL_.*
YLabel=$\mathrm{d}\Delta\sigma_{LL}/\mathrm{d}X$ [pb / unit]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/Sigma(UU|PP|PM|MP|MM)_.*
YLabel=$\mathrm{d}\sigma/\mathrm{d}X$ [pb / unit]
LogY=1
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/Shape(UU|PP|PM|MP|MM)_.*
YLabel=Unit-normalized shape
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(A2UU|A2LL)_.*
XLabel=Integrated angular moment
YLabel=$2\langle\cos(2\psi)\rangle$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(B2UU|B2LL)_.*
XLabel=Integrated angular moment
YLabel=$2\langle\sin(2\psi)\rangle$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_dpsi12_.*
XLabel=$\Delta\psi_{12}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_hardplane_primary_.*
XLabel=$\psi_{\mathrm{hard},1}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_interjet_dpsi11_.*
XLabel=$\Delta\psi_{11'}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_eeec_squeezed_.*
XLabel=$\phi_{(ij)k}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_bz_angle$
XLabel=$\chi_{\mathrm{BZ}}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_q2_beta(1|2)_.*
XLabel=$Q_{2,\beta}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_s2_beta(1|2)_.*
XLabel=$S_{2,\beta}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_(dijet_threshold_denominator|ge3_threshold|ge4_threshold)$
XLabel=Jet $p_T$ threshold [GeV]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(R32|R43|ThirdJetVeto)_(UU|PP|PM|MP|MM)$
XLabel=Jet $p_T$ threshold [GeV]
YLabel=Rate ratio or veto efficiency
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_pt(31|41)_cumulative_tail$
XLabel=$p_{T,n}/p_{T,1}$ lower threshold
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/COMPARISON/OnMinusOff_.*
YLabel=Spin on $-$ spin off
LogY=0
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/DIAGNOSTICS/SingleSpinA_.*
YLabel=$A_L^{(1)}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/DIAGNOSTICS/SingleSpinB_.*
YLabel=$A_L^{(2)}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES/DIAGNOSTICS/Parity_.*
YLabel=Parity residual
END PLOT
