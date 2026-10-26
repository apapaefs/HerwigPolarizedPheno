# BEGIN PLOT /PHENIX_2023_I2033856/d01-x01-y01
Title=PHENIX inclusive prompt-photon invariant cross section
XLabel=$p_T$ [GeV]
YLabel=$E\mathrm{d}^3\sigma/\mathrm{d}p^3$ [pb GeV$^{-2}$]
LogY=1
RatioPlot=1
RatioPlotYLabel=$(\mathrm{MC}-\mathrm{Data})/\sigma_{\mathrm{Data}}$
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/d01-x01-y02
Title=PHENIX isolated prompt-photon invariant cross section
XLabel=$p_T$ [GeV]
YLabel=$E\mathrm{d}^3\sigma/\mathrm{d}p^3$ [pb GeV$^{-2}$]
LogY=1
RatioPlot=1
RatioPlotYLabel=$(\mathrm{MC}-\mathrm{Data})/\sigma_{\mathrm{Data}}$
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/d02-x01-y01$
Title=PHENIX isolated prompt-photon double-spin asymmetry
XLabel=$p_T$ [GeV]
YLabel=$A_{LL}$
RatioPlot=1
RatioPlotYLabel=$(\mathrm{MC}-\mathrm{Data})/\sigma_{\mathrm{Data}}$
LogY=0
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/Parity_*
Title=Prompt-photon parity closure
XLabel=$p_T$ [GeV]
YLabel=normalized residual
RatioPlot=0
LogY=0
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/Pull_*
Title=PHENIX data--prediction pull
XLabel=$p_T$ [GeV]
YLabel=$(\mathrm{theory}-\mathrm{data})/\sigma$
RatioPlot=0
LogY=0
YMin=-3.5
YMax=3.5
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/DIAGNOSTICS/.*
XLabel=$p_T$ [GeV]
RatioPlot=0
LogY=0
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/DIAGNOSTICS/Component_.*_cross_section$
Title=PHENIX prompt photons: individual generated contribution
YLabel=$E\mathrm{d}^3\sigma/\mathrm{d}p^3$ [pb GeV$^{-2}$]
LogY=1
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/DIAGNOSTICS/Component_.*_isolated_all$
Title=PHENIX isolated photons: individual contribution asymmetry
YLabel=$A_{LL}$
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/DIAGNOSTICS/Parity_.*
Title=Prompt-photon parity closure
YLabel=normalized residual
# END PLOT
# BEGIN PLOT /PHENIX_2023_I2033856/DIAGNOSTICS/SingleSpin.*
Title=Prompt-photon single-spin closure
YLabel=$A_L$
# END PLOT
