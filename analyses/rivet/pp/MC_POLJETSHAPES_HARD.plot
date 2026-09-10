BEGIN PLOT /MC_POLJETSHAPES_HARD/.*
RatioPlot=0
LogY=0
ConnectBins=1
LegendAlign=l
LegendAnchor=upper left
LegendXPos=0.00
LegendYPos=0.95
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/ALL_.*
YLabel=$A_{LL}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/DeltaSigmaLL_.*
YLabel=$\mathrm{d}\Delta\sigma_{LL}/\mathrm{d}X$ [pb / unit]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/Sigma(UU|PP|PM|MP|MM)_.*
YLabel=$\mathrm{d}\sigma/\mathrm{d}X$ [pb / unit]
LogY=1
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/Shape(UU|PP|PM|MP|MM)_.*
YLabel=Unit-normalized shape
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(A2UU|A2LL)_.*
XLabel=Integrated angular moment
YLabel=$2\langle\cos(2\psi)\rangle$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(B2UU|B2LL)_.*
XLabel=Integrated angular moment
YLabel=$2\langle\sin(2\psi)\rangle$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/C2(UU|LL|PP|PM|MP|MM)_pt[0-9]+_[0-9]+_resolved_dphi31_vs_pt31$
XLabel=$p_{T,3}/p_{T,1}$
YLabel=$2\langle\cos(2\Delta\phi_{3|1})\rangle$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/S2(UU|LL|PP|PM|MP|MM)_pt[0-9]+_[0-9]+_resolved_dphi31_vs_pt31$
XLabel=$p_{T,3}/p_{T,1}$
YLabel=$2\langle\sin(2\Delta\phi_{3|1})\rangle$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/C2(UU|LL|PP|PM|MP|MM)_pt[0-9]+_[0-9]+_resolved_dpsi34_vs_pt41$
XLabel=$p_{T,4}/p_{T,1}$
YLabel=$2\langle\cos(2\Delta\psi_{34})\rangle$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/S2(UU|LL|PP|PM|MP|MM)_pt[0-9]+_[0-9]+_resolved_dpsi34_vs_pt41$
XLabel=$p_{T,4}/p_{T,1}$
YLabel=$2\langle\sin(2\Delta\psi_{34})\rangle$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_pt[0-9]+_[0-9]+_dpsi12_.*
XLabel=$\Delta\psi_{12}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_pt[0-9]+_[0-9]+_hardplane_primary_.*
XLabel=$\psi_{\mathrm{hard},1}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_pt[0-9]+_[0-9]+_interjet_dpsi11_.*
XLabel=$\Delta\psi_{11'}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_eeec_squeezed_.*
XLabel=$\phi_{(ij)k}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_bz_angle$
XLabel=$\chi_{\mathrm{BZ}}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_pt[0-9]+_[0-9]+_resolved_dphi31$
XLabel=$\Delta\phi_{3|1}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM|ShapeUU|ShapePP|ShapePM|ShapeMP|ShapeMM)_pt[0-9]+_[0-9]+_resolved_dpsi34$
XLabel=$\Delta\psi_{34}$ [rad]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet(1|2)_nconst$
XLabel=Jet constituent multiplicity $N_{\mathrm{const}}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_jet(1|2)_ptd$
XLabel=$p_T^D$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_q2_beta(1|2)_.*
XLabel=$Q_{2,\beta}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_s2_beta(1|2)_.*
XLabel=$S_{2,\beta}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_pt[0-9]+_[0-9]+_(dijet_threshold_denominator|ge3_threshold|ge4_threshold)$
XLabel=Jet $p_T$ threshold [GeV]
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(R32|R43|ThirdJetVeto)_(UU|PP|PM|MP|MM)_pt[0-9]+_[0-9]+$
XLabel=Jet $p_T$ threshold [GeV]
YLabel=Rate ratio or veto efficiency
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/(ALL|DeltaSigmaLL|SigmaUU|SigmaPP|SigmaPM|SigmaMP|SigmaMM)_pt[0-9]+_[0-9]+_pt(31|41)_cumulative_tail$
XLabel=$p_{T,n}/p_{T,1}$ lower threshold
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/COMPARISON/OnMinusOff_.*
YLabel=Full spin $-$ LHE-like
LogY=0
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/DIAGNOSTICS/SingleSpinA_.*
YLabel=$A_L^{(1)}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/DIAGNOSTICS/SingleSpinB_.*
YLabel=$A_L^{(2)}$
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/DIAGNOSTICS/Parity_.*
YLabel=Parity residual
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/.*_pt20_30_.*
Title=Particle jet $20\leq p_T<30$ GeV
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/.*_pt30_45_.*
Title=Particle jet $30\leq p_T<45$ GeV
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/.*_pt20_30_jet[34]_pt$
Title=Leading jet $20\leq p_T<30$ GeV
XMin=2
XMax=30
END PLOT

BEGIN PLOT /MC_POLJETSHAPES_HARD/.*_pt30_45_jet[34]_pt$
Title=Leading jet $30\leq p_T<45$ GeV
XMin=2
XMax=45
END PLOT
