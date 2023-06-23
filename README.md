# ConeTransferMatrix
Transfer Matrix implementations for Conical Waveguides. The code in this repository accompanies an arxiv preprint by Grothe et al.
"A Transfer Matrix for the Input Impedance of weakly tapered Cones as of Wind Instruments" (https://arxiv.org/abs/2303.12750). It can be used to produce the results presented in Fig. 2.

How to use:
Run main_compare_conemodels.m in Matlab to produce a comparison of impedance results obtained with the transmission matrix method, against a reference curve obtained with 1D finite-elements (closed_cone(ow).txt).
Required Matlab functions: cons_ow.m, TM_con_kulik.m, TM_con_nederveen.m, TM_cyl.m, TML.m
Reference data: closed_cone(ow).txt 
The reference data is computed with the python library openWind, using the script closed_cone.py (kindly provided by Augustin Ernoult).
