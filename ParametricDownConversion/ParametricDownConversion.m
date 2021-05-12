(* ::Package:: *)

(* Wolfram Language Package *)

(* Created by the Wolfram Workbench 20 Oct 2018 *)

BeginPackage["ParametricDownConversion`"]
(* Exported symbols added here with SymbolName::usage *) 

FWHMextremepoints::usage="FWHMextremepoints[input_]
Find FWHM extreme points of a list of data of the form {{x,f(x)},...}. "

FWHM::usage="FWHM[data_]
Find FWHM of a list of data of the form {{x,f(x)},...}."

MiddlePosition::usage="MiddlePosition[data_]
Find middle position respect to the FWHM of a list of data of the form {{x,f(x)},...}."

SoL::usage="SoL=299792458
Speed of light."

LambdaToOmega::usage="LambdaToOmega[l_]
Convert wavelength to angular frequency: \[Lambda]=2*\[Pi]*c/\[Omega]."

OmegaToLambda::usage="OmegaToLambda[o_]
Convert angular frequency to wavelength: \[Omega]=2*\[Pi]*c/\[Lambda]."

DeltaLambdaToDeltaOmega::usage="DeltaLambdaToDeltaOmega[dl_,l_]
Convert a wavelength range to an angular frequency range: ((c*\[CapitalDelta]\[Lambda])/\[Lambda]p^2)*2*\[Pi].
Takes as arguments \[CapitalDelta]\[Lambda] and central \[Lambda]"

DeltaOmegaToDeltaLambda::usage="DeltaOmegaToDeltaLambda[do_,l_]
Convert an angular frequency range to a wavelength range: (\[Lambda]^2*\[CapitalDelta]\[Omega])/(c*2*\[Pi]).
Takes as arguments \[CapitalDelta]\[Omega] and central \[Lambda]."

Wavevector::usage="Wavevector[l_,n_:Unitize,T___]
Find the wavevector k=(2\[Pi]/\[Lambda])*n(\[Lambda]).
Takes as arguments \[Lambda] and the refractive index function n[\[Lambda]]. Default n[\[Lambda]] value is 1. 
It takes a third optional argument(s) which is/are the refractive index parameters (usually the temperature if needed)."

DeltaK::usage="DeltaK[l1_,l2_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]
Compute the momentum mismatch in the second order nonlinear process.
Takes as arguments \[Lambda]1,\[Lambda]2,n1,n2,n3, where n1,2,3 are the refractive indeces of the three photons (signal, idler, pump).
It takes a \!\(\*SuperscriptBox[\(6\), \(th\)]\) optional argument(s) which is/are the refractive index parameters (usually the temperature if needed).
E.g. In type two down conversion use: 
DeltaK[\[Lambda]1,\[Lambda]2,KTPRefIndY,KTPRefIndZ,KTPRefIndY]"

GroupVelocityInverse::usage="GroupVelocityInverse[l1_,l2_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]
Compute the inverse group velocity of the first wavelength taken as argument.
Takes as arguments \[Lambda]1,\[Lambda]2,n1,n2,n3, where n1,2,3 are the refractive indeces of the three photons (signal, idler, pump).
It takes a \!\(\*SuperscriptBox[\(6\), \(th\)]\) optional argument(s) which is/are the refractive index parameters (usually the temperature if needed).
"

DispersionParameter::usage="DispersionParameter[l1_,l2_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]
Compute the dispersion parameter D and the angle \[Theta] in the (\!\(\*SubscriptBox[\(\[Omega]\), \(s\)]\),\!\(\*SubscriptBox[\(\[Omega]\), \(i\)]\)) plane in rads and degrees.
Takes as arguments \[Lambda]1,\[Lambda]2,n1,n2,n3, where n1,2,3 are the refractive indeces of the three photons (signal, idler, pump).
It takes a \!\(\*SuperscriptBox[\(6\), \(th\)]\) optional argument(s) which is/are the refractive index parameters (usually the temperature if needed).
"

KTPRefIndY::usage="KTPRefIndY[l_]
Refractive index KTP on Y axis."

KTPRefIndZ::usage="KTPRefIndZ[l_]
Refractive index KTP on Z axis."

KTPDeltaKtype2::usage="KTPDeltaKtype2[l1_,l2_]
Find DeltaK for the KTP type2 PDC."

KTPQPMVector::usage="KTPQPMVector[pp_, T_:40]
Compute quasi phase matching vector for the KTP.
Takes as arguments the poling period and the temperature."

KTPPolingPeriod::usage="KTPPolingPeriod[\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,T_:40,n1_:KTPRefIndY,n2_:KTPRefIndZ,n3_:KTPRefIndY]
Find the poling period for the ppKTP.
Takes as arguments \[Lambda]1,\[Lambda]2,\[Lambda]3,T,n1,n2,n3. The default value for T is 40, the value for n1,n2,n3 are the standard values for the KTP."

MgOLNRefIndE::usage="MgOLNRefIndE[\[Lambda]_,T_]
Refractive index of MgO:LN on E axis, from Covesion website."

LNRefIndE::usage="LNRefIndE[\[Lambda]_,T_]
Refractive index of LN on E axis, from Covesion website."

SpectralAMPfwhmGauss::usage="SpectralAMPfwhmGauss[\[Sigma]_]
FWHM spectral amplitude of the function Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))] given \[Sigma] parameter."

SpectralINTfwhmGauss::usage="SpectralINTfwhmGauss[\[Sigma]_]
FWHM spectral intensity of the function Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))] given \[Sigma] parameter."

SpectralAMPfwhmSech::usage="SpectralAMPfwhmSech[\[Tau]_]
FWHM spectral amplitude of the function Sech[\[Pi] x \[Tau] / 2] given \[Tau] parameter."

SpectralINTfwhmSech::usage="SpectralINTfwhmSech[\[Tau]_]
FWHM spectral intensity of the function Sech[\[Pi] x \[Tau] / 2] given \[Tau] parameter."

SpectralAMPfwhmSinc::usage="SpectralAMPfwhmSinc[\[Delta]_] 
FWHM spectral amplitude of the function Sinc[\[Delta]*x] given \[Delta] parameter."

SpectralINTfwhmSinc::usage="SpectralINTfwhmSinc[\[Delta]_]
FWHM spectral intensity of the function Sinc[\[Delta]*x] given \[Delta] parameter."

TemporalAMPfwhmGauss::usage="TemporalAMPfwhmGauss[\[Sigma]_]
FWHM temporal amplitude of the function Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))] given \[Sigma] parameter."

TemporalINTfwhmGauss::usage="TemporalINTfwhmGauss[\[Sigma]_]
FWHM temporal intensity of the function Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))] given \[Sigma] parameter."

TemporalAMPfwhmSech::usage="TemporalAMPfwhmSech[\[Tau]_]
FWHM temporal amplitude of the function Sech[\[Pi] x \[Tau] / 2] given \[Tau] parameter."

TemporalINTfwhmSech::usage="TemporalINTfwhmSech[\[Tau]_]
FWHM temporal intensity of the function Sech[\[Pi] x \[Tau] / 2] given \[Tau] parameter."

TemporalAMPfwhmSinc::usage="TemporalAMPfwhmSinc[\[Delta]_]
FWHM temporal amplitude of the function Sinc[\[Delta]*x] given \[Delta] parameter."

TemporalINTfwhmSinc::usage="TemporalINTfwhmSinc[\[Delta]_]
FWHM temporal intensity of the function Sinc[\[Delta]*x] given \[Delta] parameter."

SpectralAMPGaussToSech::usage="SpectralAMPGaussToSech[\[Sigma]_]
Given \[Sigma] find \[Tau] that matches the spectral amplitudes FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

SpectralAMPSechToGauss::usage="SpectralAMPSechToGauss[\[Tau]_]
Given \[Tau] find \[Sigma] that matches the spectral amplitudes FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

SpectralAMPGaussToSinc::usage="SpectralAMPGaussToSinc[\[Sigma]_]
Given \[Sigma] find \[Delta] that matches the spectral amplitudes FWHM of the functions Sinc[\[Delta]*x] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

SpectralAMPSincToGauss::usage="SpectralAMPSincToGauss[\[Delta]_]
Given \[Delta] find \[Sigma] that matches the spectral amplitudes FWHM of the functions Sinc[\[Delta]*x] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

SpectralAMPSechToSinc::usage="SpectralAMPSechToSinc[\[Tau]_]
Given \[Tau] find \[Delta] that matches the spectral amplitudes FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Sinc[\[Delta]*x]."

SpectralAMPSincToSech::usage="SpectralAMPSincToSech[\[Delta]_]
Given \[Delta] find \[Tau] that matches the spectral amplitudes FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Sinc[\[Delta]*x]."

SpectralINTGaussToSech::usage="SpectralINTGaussToSech[\[Sigma]_]
Given \[Sigma] find \[Tau] that matches the spectral intensities FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

SpectralINTSechToGauss::usage="SpectralINTSechToGauss[\[Tau]_]
Given \[Tau] find \[Sigma] that matches the spectral intensities FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

SpectralINTGaussToSinc::usage="SpectralINTGaussToSinc[\[Sigma]_]
Given \[Sigma] find \[Delta] that matches the spectral intensities FWHM of the functions Sinc[\[Delta]*x] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

SpectralINTSincToGauss::usage="SpectralINTSincToGauss[\[Delta]_]
Given \[Delta] find \[Tau] that matches the spectral intensities FWHM of the functions Sinc[\[Delta]*x] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

SpectralINTSincToSech::usage="SpectralINTSechToSinc[\[Delta]_]
Given \[Delta] find \[Tau] that matches the spectral intensities FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Sinc[\[Delta]*x]."

SpectralINTSechToSinc::usage="SpectralINTSincToSech[\[Delta]_]
Given \[Tau] find \[Delta] that matches the spectral intensities FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Sinc[\[Delta]*x]."

TemporalAMPGaussToSech::usage="TemporalAMPGaussToSech[\[Sigma]_]
Given \[Sigma] find \[Tau] that matches the temporal amplitudes FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

TemporalAMPSechToGauss::usage="TemporalAMPSechToGauss[\[Tau]_]
Given \[Tau] find \[Sigma] that matches the temporal amplitudes FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

TemporalAMPGaussToSinc::usage="TemporalAMPSincToGauss[\[Sigma]_]
Given \[Sigma] find \[Delta] that matches the temporal amplitudes FWHM of the functions Sinc[\[Delta]*x] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

TemporalAMPSincToGauss::usage="TemporalAMPSincToGauss[\[Delta]_]
Given \[Tau] find \[Sigma] that matches the temporal amplitudes FWHM of the functions Sinc[\[Delta]*x] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

TemporalAMPSincToSech::usage="TemporalAMPSincToSech[\[Delta]_]
Given \[Delta] find \[Tau] that matches the temporal amplitudes FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Sinc[\[Delta]*x]."

TemporalAMPSechToSinc::usage="TemporalAMPSechToSinc[\[Tau]_]
Given \[Tau] find \[Delta] that matches the temporal amplitudes FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Sinc[\[Delta]*x]."

TemporalINTGaussToSech::usage="TemporalINTGaussToSech[\[Sigma]_]
Given \[Sigma] find \[Tau] that matches the temporal intensities FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

TemporalINTSechToGauss::usage="TemporalINTSechToGauss[\[Tau]_]
Given \[Tau] find \[Sigma] that matches the temporal intensities FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

TemporalINTGaussToSinc::usage="TemporalINTGaussToSinc[\[Sigma]_]
Given \[Sigma] find \[Delta] that matches the temporal intensities FWHM of the functions Sinc[\[Delta]*x] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

TemporalINTSincToGauss::usage="TemporalINTSincToGauss[\[Delta]_]
Given \[Delta] find \[Sigma] that matches the temporal intensities FWHM of the functions Sinc[\[Delta]*x] and Exp[-\!\(\*SuperscriptBox[\(x\), \(2\)]\)/(2 \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\))]."

TemporalINTSincToSech::usage="TemporalINTSincToSech[\[Delta]_]
Given \[Delta] find \[Tau] that matches the temporal intensities FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Sinc[\[Delta]*x]."

TemporalINTSechToSinc::usage="TemporalINTSechToSinc[\[Tau]_]
Given \[Tau] find \[Delta] that matches the temporal intensities FWHM of the functions Sech[\[Pi] x \[Tau] / 2] and Sinc[\[Delta]*x]."


SigmaFromPulseLength::usage="SigmaFromPulseLength[pulseLength_]
Given pulse length (temporal intensity FWHM) find \[Sigma] parameter."

TauFromPulseLength::usage="TauFromPulseLength[pulseLength_]`
Given pulse length (temporal intensity FWHM) find \[Tau] parameter."

GaussToSechPulse::usage="GaussToSechPulse[pulseLength_]
Given pulse length (temporal intensity FWHM of Gaussian) return pulse length of \!\(\*SuperscriptBox[\(Sech\), \(2\)]\) pulse having the same spectral amplitude FWHM."

SechToGaussPulse::usage="SechToGaussPulse[pulseLength_]
Given pulse length (temporal intensity FWHM of \!\(\*SuperscriptBox[\(Sech\), \(2\)]\)) return pulse length of Gaussian pulse having the same spectral amplitude FWHM."


PEFGauss::usage="PEFGauss[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Sigma]_]
Gaussian pump envelope function: Exp[-(((\[Omega]1 + \[Omega]2) - \[Omega]0)^2/(2*\[Sigma]^2))]
Takes as arguments \[Omega]1 and \[Omega]2, \[Omega]0 and \[Sigma]."

PEFSech::usage="PEFSech[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Tau]_]
Sech pump envelope function: Sech[(Pi*\[Tau]*((\[Omega]1 + \[Omega]2) - \[Omega]0))/2]
Takes as arguments \[Omega]1 and \[Omega]2, \[Omega]0 and \[Tau]."

PEFGaussChirp::usage="PEFGaussChirp[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Sigma]_, k_]
Gaussian pump envelope function with chirping: Exp[-(((\[Omega]1 + \[Omega]2) - \[Omega]0)^2/(2*\[Sigma]^2))]Exp[-\[ImaginaryI]*k*((\[Omega]1 + \[Omega]2) - \[Omega]0)^2]
Takes as arguments \[Omega]1 and \[Omega]2, \[Omega]0, \[Sigma] and k, where k is the chirping parameter."

PEFSechChirp::usage="PEFSechChirp[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Tau]_, k_]
Sech pump envelope function with chirping: Sech[(Pi*\[Tau]*((\[Omega]1 + \[Omega]2) - \[Omega]0))/2]Exp[-\[ImaginaryI]*k*((\[Omega]1 + \[Omega]2) - \[Omega]0)^2]
Takes as arguments \[Omega]1 and \[Omega]2, \[Omega]0, \[Tau] and k, where k is the chirping parameter."

PMFSinc::usage="PEFGauss[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Delta]_]
Sinc PMF. Takes as arguments \[Omega]1 and \[Omega]2, \[Omega]0 and \[Delta]."

PMFGauss::usage="PMFGauss[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Sigma]_]
Gaussian PMF. Takes as arguments \[Omega]1 and \[Omega]2, \[Omega]0 and \[Sigma]."

PMFAntisymmetric::usage="PMFAntisymmetric[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Sigma]_]
Antisymmetric PMF. Takes as arguments \[Omega]1 and \[Omega]2, \[Omega]0 and \[Sigma]."


FrequencyList::usage="FrequencyList[min_,max_,discretisation_]
Takes as argument the min,man and step values of the frequencies, and creates a list."

FrequencyMatrix::usage="FrequencyMatrix[frequencies1_,frequencies2_]
Takes two list of frequencies and combine them."

MatrixToLambdas::usage="MatrixToLambdas[matrix_,scaling___]
Convert a list of the form {{\[Omega]a,\[Omega]b,f(\[Omega]a,\[Omega]b)},...} to frequencies {{\[Lambda]a,\[Lambda]b,f(\[Omega]a,\[Omega]b)},...}.
Scaling can be applied (as a multiplicative factor) to the lambdas."

PEFMatrix::usage="PEFMatrix[frequencyMatrix_,pumpFunction_,param__]
Calculate the PEF, having as input the list of signal and idler combinations.
It can also be used to calculate the PMF if we have it's analytical form instead of the list of domains."

CrystalWalls::usage = "CrystalWalls[widths_]
crystalWalls[widths] takes a list of domain widths as input and finds the domain walls in a crystal.
The output is of the form: {0,w1,w2,w3,...,wn}"

ImportCrystal::usage="ImportCrystal[file_, directory_:NotebookDirectory[]]
Import crystal from input file where the domain widths are given, and finds walls position. Takes as argument file name and directory name."

ImportReverseCrystal::usage="ImportReverseCrystal[file_, directory_:NotebookDirectory[]]
Import reverse crystal from input file where the domain widths are given, and finds walls position. Takes as argument file name and directory name."

PeriodicallyPoledEquivalentPMF::usage="PeriodicallyPoledEquivalentPMF[dk_, polingperiod_, crystallength_]
Calculate the equivalent PMF of a periodically poled crystal at a momentum mismatch dk given the poling period and the crystal length."

ApodizedGaussianEquivalentPMF::usage="ApodizedGaussianEquivalentPMF[dk_, polingperiod_, crystallength_, \[Sigma]den_]
Calculate the equivalent PMF of an apodized Gaussian crystal at a momentum mismatch dk given the poling period, the crystal length and the width parameter."

PeriodicallyPoledCrystal::usage="PeriodicallyPoledCrystal[polingperiod_, crystallength_]
Build periodically poled crystal having as input poling period and crystal length."

ApodizedGaussianCrystal::usage="ApodizedGaussianCrystal[polingperiod_, subdivision_, crystallength_, \[Sigma]den_, exportOrientation_:False]
Build apodized Gaussian poled crystal having as input poling period, number of subdivisions of the standard domain, crystal length and sigma parameter.
If the last optional argument is True exports also the domain orientation structure of the crystal."

ApodizedAntisymmetricCrystal::usage="ApodizedAntisymmetricCrystal[polingperiod_, subdivision_, crystallength_, \[Sigma]den_, exportOrientation_:False]
Build apodized antisymmetric poled crystal having as input poling period, number of subdivisions of the standard domain, crystal length and sigma parameter.
If the last optional argument is True exports also the domain orientation structure of the crystal."

ApodizedGenericCrystal::usage="ApodizedGenericCrystal[polingperiod_, subdivision_, crystallength_, target_,  exportOrientation_:False]
Build apodized generic poled crystal having as input poling period, number of subdivisions of the standard domain, crystal length and PMF at \!\(\*SubscriptBox[\(\[CapitalDelta]k\), \(0\)]\) to be tracked.
If the last optional argument is True exports also the domain orientation structure of the crystal."

DomainOrientationToDomainWidths::usage="DomainOrientationToDomainWidths[orientation_,domWidth_]
Build domain widths having as input orientation of each domain and domain width."

Amplitude::usage="Amplitude[z_, \[CapitalDelta]k_, domainList_]
Return amplitude at position z in the crystal and DeltaK value, takes as input position z, DeltaK value and domain wall positions"

PhaseMatchingFunction::usage="PhaseMatchingFunction[\[CapitalDelta]k_, domainList_]
Compute PMF at given DeltaK value. Takes as argument DeltaKand the domain wall positions."

PMFBandwidth::usage="PMFBandwidth[frequencies1_,frequencies2_,domainList_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]
Calculate the PMF bandwidth (intensity) in Hz, given the frequency lists, the domain walls positions and the Sellmeiers of the material.
Function still under testing."

DeltaKMatrix::usage="DeltaKMatrix[frequencyMatrix_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]
Compute the DeltaK matrix for the photon freqencies. Takes as argument  the frequency matrix and the refractive indeces of the three photons."

PMFMatrix::usage="Compute the PMF matrix, takes as argument a DeltaK Matrix and the list of domain wall positions."

JSAMatrix::usage="JSAMatrix[PEF_, PMF_]
Compute JSA matrix given PEF and PMF matrix"

JSIMatrix::usage="JSIMatrix[JSA_]
Compute JSi matrix given JSA matrix."

JSAAbsMatrixFromJSI::usage="JSAAbsMatrixFromJSI[JSI_]
Compute absoute value of JSA matrix given JSI matrix"

JSAAbsMatrix::usage="JSAAbsMatrix[JSA_]
Compute absoute value of JSA matrix given JSA matrix"

Brightness::usage="Brightness[JSA_]
Compute brightness given JSA matrix."

DisassembleMatrix::usage="DisassembleMatrix[matrix_]
Takes as argument a list of the form {{coordA1,coordB1,value11},...} and returns three lists {{coords1},{coords2},{matrix_values}} "

ReassembleMatrix::usage="ReassembleMatrix[frequencyMatrix_,matrix_]
"

SchmidtDecomposition::usage="SchmidtDecomposition[matrix_]
Perform Schimdt Decomposition on a matrix of the form {{coordA1,coordB1,value11}.
It returns a list {purity, schmidt coefficients, signalModes,idlerModes}"

MarginalSpectra::usage="MarginalSpectra[JSA_]
Compute marginal spectrum of {signal,idler} given JSA or JSI."

ComputeFilter::usage="ComputeFilter[frequencyMatrix_,filteringFunction_,param___]
Filtering function, it returns the filter applied to signal and idler photon.
Takes frequency matrix, filtering function and optinal parameters as argument."

ComputeFilteredJSA::usage="ComputeFilteredJSA[JSA_, filterSignal_, filterIdler_]
Function for computing the filtered JSA given the JSA and the filters; in the density plot the signal frequencies on the x axis, idler frequencies on the y axis."

ComputeHeraldingEfficiency::usage="ComputeHeraldingEfficiency[JSA_, filterSignal_, filterIdler_]
Find symmetric heralding effinciency when filtering is applied to both the photons."

OptimisePulse::usage="OptimisePulse[frequencyMatrix_, PMF_, start_, stop_, step_, pumpFunction_, param__ ]
Optimise the pulse length respect to the crystal PMF."

OptimisePulseFiltered::usage="OptimisePulseFiltered[frequencyMatrix_, PMF_,filterSignal_, filterIdler_, start_, stop_, step_, pumpFunction_, param__ ]
Optimise the pulse length respect to the crystal PMF when filtering is applied."

SignalIdlerHOMDip::usage="SignalIdlerHOMDip[t_, decomposedJSA_]
Return HOM dip value at time t (accepts also lists), it only works for square matrices with the same spectral ranges for the whole dip run."

HeraldedHOMDip::usage="HeraldedHOMDip[t_, decomposedJSA1_, index1_, decomposedJSA2_, index2_]
HOM dip value at time t (accepts also lists) between two heralded photons produced by PDC process, index 1 and 2 can be either 3 or \"s\" (signal) and 4 \"i\" (idler) depending on the photon that are interfered"


packageDirectory=DirectoryName[$InputFileName];

TutorialNotebook[] :=
    SetOptions[
      NotebookOpen@FileNameJoin[{packageDirectory, "Tutorial.nb"}],
      Saveable -> False
    ]
    
Print[
"Parametric Down Conversion is a package (mainly) for simulating second order nonlinear processes, in particular PDC. 
A list of the functions in the library is shown below. 
More details on how to use the library can be found in the tutorial notebook by calling the function: TutorialNotebook[]."
]
Print["A list of of the functions can be found by calling the command: ?\"ParametricDownConversion`*\""]


Begin["`Private`"]
(* Implementation of the package *)


FWHMextremepoints[input_] := Block[{data, clip, min, max, x, y}, 
     data = input; 
     data[[All,2]] = data[[All,2]]/Max[data[[All,2]]]; 
     clip = SparseArray[Clip[data[[All,2]], {0.5, 0.5}, {0, 1}]]["ColumnIndices"]; 
     If[Length[clip] < 1, Return[$Failed]]; 
     {min, max} = clip[[{1, -1},1]]; 
     min = y /. First[Solve[data[[min - 1]] + x*(data[[min]] - data[[min - 1]]) == {y, 0.5}, {x, y}]]; 
     max = y /. First[Solve[data[[max]] + x*(data[[max + 1]] - data[[max]]) == {y, 0.5}, {x, y}]]; 
     {min, max}]
 
FWHM[data_] := Block[{temp, fwhm, peak}, 
	temp = FWHMextremepoints[data]; 
     Abs[temp[[-1]] - temp[[1]]]]
 
MiddlePosition[data_] := Block[{temp, fwhm, peak}, 
	temp = FWHMextremepoints[data]; 
	Mean[temp]]


SoL=299792458;


LambdaToOmega[l_]:=2*Pi*(SoL/l);
OmegaToLambda[o_]:=2*Pi*(SoL/o);

DeltaLambdaToDeltaOmega[dl_,l_]:=((SoL*dl)/l^2)*2*Pi;
DeltaOmegaToDeltaLambda[do_,l_]:=(l^2*do)/(SoL*2*Pi);


Wavevector[l_,n_:Unitize,T___]:=(2*(Pi/l))*n[l,T];

DeltaK[l1_,l2_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]:=
	Wavevector[(1/l1 + 1/l2)^(-1),n3,T]-Wavevector[l1,n1,T]-Wavevector[l2,n2,T]
	
GroupVelocityInverse[l1_,l2_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]:=Block[
	{omega,k3,k1,k2,dk,gvi},
	k3[omega_]:=(omega+LambdaToOmega[l2])/SoL*n3[OmegaToLambda[omega+LambdaToOmega[l2]],T];
	k1[omega_]:=omega/SoL*n1[OmegaToLambda[omega],T];
	k2=(2*Pi)/l2*n2[l2,T];
	dk[omega_]:=k3[omega]-k1[omega]-k2;
	gvi[omega_]=D[dk[omega],omega];
	gvi[OmegaToLambda[l1]]	
	]
	
DispersionParameter[l1_,l2_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]:=Block[
	{gvi1,gvi2},
	gvi1 = GroupVelocityInverse[l1,l2,n1,n2,n3,T];
	gvi2 = GroupVelocityInverse[l2,l1,n2,n1,n3,T];
	{-gvi1/gvi2,-ArcTan[gvi1/gvi2],-ArcTan[gvi1/gvi2]*180/\[Pi]}
	]


KTPRefIndY[l_]:=Block[{Ay1 = 2.0993,Ay2 = 0.922683,Ay3 = 0.0467695,Ay4 = 0.0138408},
	Sqrt[Ay1 + Ay2/(1 - Ay3/(l*10^6)^2) - Ay4*(l*10^6)^2]];

KTPRefIndZ[l_]:=Block[{Az = 2.12725,Bz = 1.18431,Cz = 0.05148520000000001,Dz = 0.6603,Ez = 100.00507,Fz = 0.00968956},
	Sqrt[Az + Bz/(1 - Cz/(l*10^6)^2)+Dz/(1 - Ez/(l*10^6)^2) - Fz*(l*10^6)^2]]

KTPDeltaKtype2[l1_,l2_]:=DeltaK[l1,l2,KTPRefIndY,KTPRefIndZ,KTPRefIndY]

KTPQPMVector[pp_, T_:40] := Block[{\[Alpha]d = 6.7*^-6,\[Beta]d = 11*10^-9},
	(2*Pi)/(pp*(1 + \[Alpha]d*(T - 25) + \[Beta]d*(T - 25)^2))]

KTPPolingPeriod[\[Lambda]1_:1550*10^-9,\[Lambda]2_:1550*10^-9,\[Lambda]3_:775*10^-9,T_:40,n1_:KTPRefIndY,n2_:KTPRefIndZ,n3_:KTPRefIndY] := 
    Abs[NSolve[{Wavevector[\[Lambda]3,n3] - Wavevector[\[Lambda]1,n1] - Wavevector[\[Lambda]2,n2] - KTPQPMVector[x, T] == 0}, x][[1,1,2]]]


MgOLNRefIndE[\[Lambda]_,T_]:=Block[
	{f,a1=5.756,a2=0.0983,a3=0.2020,a4=189.32,a5=12.52,a6=1.32*10^-2,b1=2.860*10^-6,b2=4.700*10^-8,b3=6.113*10^-8,b4=1.516*10^-4},
	f=(T-24.5)*(T+570.82);
	Sqrt[a1+b1*f+(a2+b2*f)/(\[Lambda]^2-(a3+b3*f)^2)+(a4+b4*f)/(\[Lambda]^2-a5^2)-a6*\[Lambda]^2]
]

LNRefIndE[\[Lambda]_,T_]:=Block[
	{f,a1=5.35583,a2=0.100473,a3=0.20692,a4=100,a5=11.34927,a6=1.5334*10^-2,b1=4.629*10^-7,b2=3.862*10^-8,b3=-8.9*10^-9,b4=2.657*10^-5},
	f=(T-24.5)*(273.15+T+570.82);
	Sqrt[a1+b1*f+(a2+b2*f)/(\[Lambda]^2-(a3+b3*f)^2)+(a4+b4*f)/(\[Lambda]^2-a5^2)-a6*\[Lambda]^2]
]


SpectralAMPfwhmGauss[\[Sigma]_] := 2*Sqrt[2*Log[2]]*\[Sigma]
 
SpectralINTfwhmGauss[\[Sigma]_] := 2*Sqrt[Log[2]]*\[Sigma]
 
SpectralAMPfwhmSech[\[Tau]_] := ((4*ArcCosh[2])/Pi)*(1/\[Tau])
 
SpectralINTfwhmSech[\[Tau]_] := ((4*ArcSech[1/Sqrt[2]])/Pi)*(1/\[Tau])

SpectralAMPfwhmSinc[\[Delta]_] := 2*1.895494267033981`/\[Delta]
 
SpectralINTfwhmSinc[\[Delta]_] := 2*1.39155737825151`/\[Delta]

TemporalAMPfwhmGauss[\[Sigma]_] := 2*Sqrt[2*Log[2]]*(1/\[Sigma])
 
TemporalINTfwhmGauss[\[Sigma]_] := 2*Sqrt[Log[2]]*(1/\[Sigma])
 
TemporalAMPfwhmSech[\[Tau]_] := 2*ArcSech[1/2]*\[Tau]
 
TemporalINTfwhmSech[\[Tau]_] := 2*ArcSech[1/Sqrt[2]]*\[Tau]

TemporalAMPfwhmSinc[\[Delta]_] := 2*\[Delta]
 
TemporalINTfwhmSinc[\[Delta]_] := 2*\[Delta]
 
SpectralAMPGaussToSech[\[Sigma]_] := (ArcCosh[2]*Sqrt[2/Log[2]])/(Pi*\[Sigma])

SpectralAMPSechToGauss[\[Tau]_] := (ArcCosh[2]*Sqrt[2/Log[2]])/(Pi*\[Tau])

SpectralAMPGaussToSinc[\[Sigma]_] := 1.60988460331292`/\[Sigma]

SpectralAMPSincToGauss[\[Delta]_] := 1.60988460331292`/\[Delta]

SpectralAMPSechToSinc[\[Tau]_] := 2.260843295803358` \[Tau]

SpectralAMPSincToSech[\[Delta]_] := 0.4423128316129775` \[Delta]
 
SpectralINTGaussToSech[\[Sigma]_] := (2*ArcSech[1/Sqrt[2]])/(Pi*\[Sigma]*Sqrt[Log[2]])

SpectralINTSechToGauss[\[Tau]_] := (2*ArcSech[1/Sqrt[2]])/(Pi*\[Tau]*Sqrt[Log[2]])

SpectralINTGaussToSinc[\[Sigma]_] := 1.6714307501300105`/\[Sigma]

SpectralINTSincToGauss[\[Delta]_] := 1.6714307501300105`/\[Delta]

SpectralINTSechToSinc[\[Tau]_] := 2.4800530109751753` \[Tau]

SpectralINTSincToSech[\[Delta]_] := 0.4032171875256782` \[Delta]

TemporalAMPGaussToSech[\[Sigma]_] := Sqrt[2*Log[2]]/(\[Sigma]*ArcSech[1/2])

TemporalAMPSechToGauss[\[Tau]_] := Sqrt[2*Log[2]]/(\[Tau]*ArcSech[1/2])

TemporalAMPGaussToSinc[\[Sigma]_] := 1.1774100225154747`/\[Sigma]

TemporalAMPSincToGauss[\[Delta]_] := 1.1774100225154747`/\[Delta]

TemporalAMPSechToSinc[\[Tau]_] := 1.3169578969248166` \[Tau]

TemporalAMPSincToSech[\[Delta]_] := 0.759325717500207` \[Delta]

TemporalINTGaussToSech[\[Sigma]_] := Sqrt[Log[2]]/(\[Sigma]*ArcSech[1/Sqrt[2]])

TemporalINTSechToGauss[\[Tau]_] := Sqrt[Log[2]]/(\[Tau]*ArcSech[1/Sqrt[2]])

TemporalINTGaussToSinc[\[Sigma]_] := 0.8325546111576977`/\[Sigma]

TemporalINTSincToGauss[\[Delta]_] := 0.8325546111576977`/\[Delta]

TemporalINTSechToSinc[\[Tau]_] := 0.8813735870195432` \[Tau]

TemporalINTSincToSech[\[Delta]_] := 1.1345926571065108` \[Delta]


SigmaFromPulseLength[pulseLength_] := (2*Sqrt[Log[2]])/pulseLength
 
TauFromPulseLength[pulseLength_] := pulseLength/(2*ArcSech[1/Sqrt[2]])

GaussToSechPulse[pulseLength_] := Block[{sigma, tau}, 
     sigma = SigmaFromPulseLength[pulseLength]; 
     tau = SpectralAMPGaussToSech[sigma]; 
     TemporalINTfwhmSech[tau]]
 
SechToGaussPulse[pulseLength_] := Block[{sigma, tau}, 
     tau = TauFromPulseLength[pulseLength]; 
     sigma = SpectralAMPSechToGauss[tau]; 
     TemporalINTfwhmGauss[sigma]]


PEFGauss[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Sigma]_] := Exp[-(((\[Omega]1 + \[Omega]2) - \[Omega]0)^2/(2*\[Sigma]^2))]
 
PEFSech[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Tau]_] := Sech[(Pi*\[Tau]*((\[Omega]1 + \[Omega]2) - \[Omega]0))/2]

PEFGaussChirp[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Sigma]_, k_] := PEFGauss[\[Omega]1, \[Omega]2, \[Omega]0, \[Sigma]]*Exp[-I*k*((\[Omega]1 + \[Omega]2) - \[Omega]0)^2]

PEFSechChirp[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Tau]_, k_] := PEFSech[\[Omega]1, \[Omega]2, \[Omega]0, \[Tau]]*Exp[-I*k*((\[Omega]1 + \[Omega]2) - \[Omega]0)^2]

PMFSinc[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Delta]_] := Sinc[((\[Omega]1 - \[Omega]2) - \[Omega]0) \[Delta]]

PMFGauss[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Sigma]_] := Exp[-(((\[Omega]1 - \[Omega]2) - \[Omega]0)^2/(2*\[Sigma]^2))]

PMFAntisymmetric[\[Omega]1_, \[Omega]2_, \[Omega]0_, \[Sigma]_] := Exp[-(((\[Omega]1 - \[Omega]2) - \[Omega]0)^2/(2*\[Sigma]^2))]*((\[Omega]1 - \[Omega]2) - \[Omega]0)


FrequencyList[min_,max_,discretisation_] := Range[min,max,discretisation]
FrequencyMatrix[frequencies1_,frequencies2_]:= Tuples[{frequencies1, frequencies2}]; 
MatrixToLambdas[matrix_,scaling___] := {Times[OmegaToLambda@#[[1]],scaling],Times[OmegaToLambda@#[[2]],scaling], #[[3]]} &/@ matrix;


PEFMatrix[frequencyMatrix_,pumpFunction_,param__]:=
	Apply[{#1, #2, pumpFunction[#1, #2, param]} & , frequencyMatrix,  {1}]


CrystalWalls[widths_] := Block[{boundaries}, 
     boundaries = Accumulate[widths]; PrependTo[boundaries, 0]]
     
ImportCrystal[file_, directory_:NotebookDirectory[]] := Block[{widths}, 
     widths = Flatten[Import[file, Path -> directory]]; CrystalWalls[widths]]
     
ImportReverseCrystal[file_, directory_:NotebookDirectory[]] := Block[{widths}, 
     widths = Flatten[Import[file, Path -> directory]]; 
     widths = Reverse[widths]; 
     CrystalWalls[widths]]


PeriodicallyPoledEquivalentPMF[dk_, polingperiod_, crystallength_] := 
	Block[{dk0=2*\[Pi]/polingperiod},
	Sinc[crystallength (dk-dk0)/2]
	]
PeriodicallyPoledEquivalentPMF[dk_List, polingperiod_, crystallength_] := {#, PeriodicallyPoledEquivalentPMF[#, polingperiod, crystallength]} & /@ dk	
	
ApodizedGaussianEquivalentPMF[dk_, polingperiod_, crystallength_, \[Sigma]den_] := 
	Block[{dk0=2*\[Pi]/polingperiod},
	(E^(-(((dk-dk0)^2 crystallength^2)/(2 \[Sigma]den^2))) Sqrt[\[Pi]/2] Sqrt[crystallength^2/\[Sigma]den^2] \[Sigma]den (Erf[(I (dk-dk0) crystallength+\[Sigma]den^2)/(Sqrt[2] \[Sigma]den)]+Erf[(-2 I (dk-dk0) crystallength+\[Sigma]den^2)/(2 Sqrt[2] \[Sigma]den)]))/crystallength
	]
ApodizedGaussianEquivalentPMF[dk_List, polingperiod_, crystallength_, \[Sigma]den_] := {#, ApodizedGaussianEquivalentPMF[#, polingperiod, crystallength, \[Sigma]den]} & /@ dk


PeriodicallyPoledCrystal[polingperiod_, crystallength_] := 
    Block[{domains}, domains = Round[2*(crystallength/polingperiod)]; 
    ConstantArray[polingperiod/2, domains]]
      
ApodizedGaussianCrystal[polingperiod_, subdivision_, crystallength_, \[Sigma]den_, exportOrientation_:False] := 
	Block[{w, domains, stddomains, length, \[Sigma], target, Amplitude, error, orientation, record, orienUP, orienDOWN, 
      errUP, errDOWN, data, flag}, 
	stddomains = Round[2*(crystallength/polingperiod)]; 
    w = polingperiod/(2*subdivision); domains = Round[crystallength/w]; 
    length = domains*w; \[Sigma] = length/\[Sigma]den; 
    target[z_] := Sqrt[2/Pi]*\[Sigma]*(Erf[length/(2*Sqrt[2]*\[Sigma])] - Erf[(length - 2*z)/(2*Sqrt[2]*\[Sigma])]); 
    Amplitude[m_, orien_] := 
       (polingperiod/(2*Pi))*(Exp[(-I)*((2*Pi)/polingperiod)*w] - 1)*
        Sum[orien[[n]]*Exp[I*((2*Pi)/polingperiod)*n*w], {n, 1, m}]; 
	error[m_, orien_] := Abs[Amplitude[m, orien] - target[m*w]]; 
    orientation = {}; 
    Table[orienUP = Append[orientation, 1]; 
          orienDOWN = Append[orientation, -1]; 
          errUP = error[m, orienUP]; 
          errDOWN = error[m, orienDOWN]; 
          If[errUP <= errDOWN, orientation = orienUP, orientation = orienDOWN]; 
     , {m, 1, domains}]; 
     data = {}; 
     flag = 1; 
     AppendTo[orientation, 0]; 
     Table[If[orientation[[i]] == orientation[[i - 1]], 
           flag++, AppendTo[data, {flag*w}]; flag = 1; ]
     , {i, 2, Length[orientation]}]; 
     If[exportOrientation,{Flatten@data,orientation[[;;-2]]},Flatten@data]]
     
ApodizedAntisymmetricCrystal[polingperiod_, subdivision_, crystallength_, \[Sigma]den_ , exportOrientation_:False] := 
	Block[{w, domains, stddomains, length, \[Sigma], target, Amplitude, error, orientation, record, orienUP, orienDOWN, 
      errUP, errDOWN, data, flag}, 
	stddomains = Round[2*(crystallength/polingperiod)]; 
    w = polingperiod/(2*subdivision); domains = Round[crystallength/w]; 
    length = domains*w; \[Sigma] = length/\[Sigma]den; 
    target[z_] := (2 E^(1/2-(length-2 z)^2/(8 \[Sigma]^2)) (-1+E^((z (-length+z))/(2 \[Sigma]^2))) \[Sigma])/\[Pi];
    Amplitude[m_, orien_] := 
       (polingperiod/(2*Pi))*(Exp[(-I)*((2*Pi)/polingperiod)*w] - 1)*
        Sum[orien[[n]]*Exp[I*((2*Pi)/polingperiod)*n*w], {n, 1, m}]; 
	error[m_, orien_] := Abs[Amplitude[m, orien] - target[m*w]]; 
    orientation = {}; 
    Table[orienUP = Append[orientation, 1]; 
          orienDOWN = Append[orientation, -1]; 
          errUP = error[m, orienUP]; 
          errDOWN = error[m, orienDOWN]; 
          If[errUP <= errDOWN, orientation = orienUP, orientation = orienDOWN]; 
     , {m, 1, domains}]; 
     data = {}; 
     flag = 1; 
     AppendTo[orientation, 0]; 
     Table[If[orientation[[i]] == orientation[[i - 1]], 
           flag++, AppendTo[data, {flag*w}]; flag = 1; ]
     , {i, 2, Length[orientation]}]; 
     If[exportOrientation,{Flatten@data,orientation[[;;-2]]},Flatten@data]]
     
     
ApodizedGenericCrystal[polingperiod_, subdivision_, crystallength_, target_,  exportOrientation_:False] := 
	Block[{w, domains, stddomains, length, \[Sigma], Amplitude, error, orientation, record, orienUP, orienDOWN, 
      errUP, errDOWN, data, flag}, 
	stddomains = Round[2*(crystallength/polingperiod)]; 
    w = polingperiod/(2*subdivision); domains = Round[crystallength/w]; 
    length = domains*w; 
    Amplitude[m_, orien_] := 
       (polingperiod/(2*Pi))*(Exp[(-I)*((2*Pi)/polingperiod)*w] - 1)*
        Sum[orien[[n]]*Exp[I*((2*Pi)/polingperiod)*n*w], {n, 1, m}]; 
	error[m_, orien_] := Abs[Amplitude[m, orien] - target[m*w]]; 
    orientation = {}; 
    Table[orienUP = Append[orientation, 1]; 
          orienDOWN = Append[orientation, -1]; 
          errUP = error[m, orienUP]; 
          errDOWN = error[m, orienDOWN]; 
          If[errUP <= errDOWN, orientation = orienUP, orientation = orienDOWN]; 
     , {m, 1, domains}]; 
     data = {}; 
     flag = 1; 
     AppendTo[orientation, 0]; 
     Table[If[orientation[[i]] == orientation[[i - 1]], 
           flag++, AppendTo[data, {flag*w}]; flag = 1; ]
     , {i, 2, Length[orientation]}]; 
     If[exportOrientation,{Flatten@data,orientation[[;;-2]]},Flatten@data]]


DomainOrientationToDomainWidths[orientation_,domWidth_]:=Block[
{fullOrientation,flag=1,data={}},
fullOrientation=Append[orientation, 0];
Table[
	If[fullOrientation[[i]] == fullOrientation[[i - 1]], 
           flag++, AppendTo[data, {flag*domWidth}]; flag = 1; ]
     , {i, 2, Length[fullOrientation]}
     ]; 
     Flatten@data
]


Amplitude[z_, \[CapitalDelta]k_, domainList_] := 
    Block[{smaller, currentDomainIndex, lastWallPosition}, 
    smaller = Flatten[Select[domainList, #1 <= z & ]]; 
    currentDomainIndex = Length[smaller]; 
    lastWallPosition = smaller[[-1]]; 
    (1/\[CapitalDelta]k)*(If[currentDomainIndex > 1,
             Sum[(-1)^(i + 1)*(Exp[I*\[CapitalDelta]k*domainList[[i]]] - Exp[I*\[CapitalDelta]k*domainList[[i - 1]]]),{i, 2, currentDomainIndex}], Unevaluated[Sequence[]]] 
             + (-1)^currentDomainIndex*(Exp[I*\[CapitalDelta]k*z] - Exp[I*\[CapitalDelta]k*lastWallPosition]))
    ]
Amplitude[z_List, \[CapitalDelta]k_List, domainList_] := Flatten[Outer[{#1, #2, Amplitude[#1, #2, domainList]} &, z, \[CapitalDelta]k], 1];

PhaseMatchingFunction[\[CapitalDelta]k_, domainList_] := 
    Block[{}, 
    (1/\[CapitalDelta]k)*Sum[(-1)^(i + 1)*(Exp[I*\[CapitalDelta]k*domainList[[i]]] - Exp[I*\[CapitalDelta]k*domainList[[i - 1]]]), 
    {i, 2, Length[domainList]}]
    ]
PhaseMatchingFunction[\[CapitalDelta]k_List, domainList_] := {#, PhaseMatchingFunction[#, domainList]} & /@ \[CapitalDelta]k

PMFBandwidth[frequencies1_,frequencies2_,domainList_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]:=
	Block[{varyingFreqs},
	varyingFreqs=Transpose[{frequencies1,frequencies2}];
	10^12*FWHM[{#[[1]]/(2\[Pi])*10^-12,Abs[PhaseMatchingFunction[DeltaK[OmegaToLambda[#[[1]]],OmegaToLambda[#[[2]]],n1,n2,n3,T],domainList]]^2}&/@varyingFreqs]
	]


DeltaKMatrix[frequencyMatrix_,n1_:Unitize,n2_:Unitize,n3_:Unitize,T___]:=
	Apply[{#1, #2,DeltaK[OmegaToLambda[#1],OmegaToLambda[#2],n1,n2,n3,T]} & , frequencyMatrix, {1}]

PMFMatrix[deltaKMatrix_,domainList_]:=
	Apply[{#1, #2, PhaseMatchingFunction[#3, domainList]} & , deltaKMatrix, {1}]


JSAMatrix[PEF_, PMF_] := If[PEF[[1 ;; All,{1, 2}]] == PMF[[1 ;; All,{1, 2}]], 
     Transpose[{PEF[[;;,1]],PEF[[;;,2]],PEF[[;;,3]]*PMF[[;;,3]]}], 
     Print["Error: Mismatch in the PMF and PEF frequencies"]]
     
JSIMatrix[JSA_] := Block[{JSI}, 
	JSI = JSA; 
	JSI[[;;,3]] = Abs[JSA[[;;,3]]]^2; 
	JSI]

JSAAbsMatrixFromJSI[JSI_] := Block[{JSA}, 
	JSA = JSI; 
	JSA[[;;,3]] = Sqrt[JSI[[;;,3]]]; 
	JSA]

JSAAbsMatrix[JSA_] := Block[{JSI}, 
	JSI = JSA; 
	JSI[[;;,3]] = Abs[JSA[[;;,3]]]; 
	JSI]
                     
Brightness[JSA_] := Total[Abs[JSA[[;;,3]]]^2]


DisassembleMatrix[matrix_]:=Block[{f1,f2,mat},
	{f1,f2}=DeleteDuplicates/@{matrix[[;;,1]],matrix[[;;,2]]};
	{f1,f2,Transpose@ArrayReshape[matrix[[;;,3]],{Length[f1],Length[f2]}]}
	]
	
ReassembleMatrix[frequencyMatrix_,matrix_]:=Block[{mat},
	mat=Flatten[Transpose[matrix]];
	Flatten[{#[[1]],#[[2]]}]&/@Transpose[{frequencyMatrix,mat}]
	]


SchmidtDecomposition[matrix_] := 
    Block[{disMatrix,decoMatrix, JSAmatrix, schmidtCoeff, purity, signalModes, 
      idlerModes}, 
      disMatrix = DisassembleMatrix[matrix]; 
      JSAmatrix = disMatrix[[3]];
      decoMatrix = SingularValueDecomposition[JSAmatrix]; 
      schmidtCoeff = Diagonal[decoMatrix[[2]]/Sqrt[Tr[decoMatrix[[2]] . decoMatrix[[2]]]]]; 
      purity = Total[schmidtCoeff^4]; 
      signalModes = (Transpose[{disMatrix[[1]], decoMatrix[[3,1 ;; All,#1]]}] & ) /@ 
        Range[Length[decoMatrix[[3]]]]; 
      idlerModes = (Transpose[{disMatrix[[2]], decoMatrix[[1,1 ;; All,#1]]}] & ) /@ 
        Range[Length[decoMatrix[[1]]]]; 
      {purity, schmidtCoeff, signalModes,idlerModes}]


MarginalSpectra[JSA_] := Block[{tempSign, tempIdl, signalSpectrum,idlerSpectrum}, 
	tempSign = Transpose[GatherBy[JSA, First], {1, 2, 3}]; 
    tempIdl = Transpose[GatherBy[JSA, First], {2, 1, 3}]; 
    signalSpectrum = ({#1[[1]]/Length[tempSign], #1[[3]]} & ) /@ Total /@ tempSign; 
    idlerSpectrum = ({#1[[2]]/Length[tempIdl], #1[[3]]} & ) /@ Total /@ tempIdl; 
    {signalSpectrum, idlerSpectrum}]


ComputeFilter[frequencyMatrix_,filteringFunction_,param___] := 
    Block[{signalFiltered, idlerFiltered}, 
     signalFiltered = Apply[{#1, #2, filteringFunction[#1,param]} & , frequencyMatrix, {1}]; 
     idlerFiltered = Apply[{#1, #2, filteringFunction[#2,param]} & , frequencyMatrix, {1}]; 
     {signalFiltered, idlerFiltered}]

ComputeFilteredJSA[JSA_, filterSignal_, filterIdler_] := 
    If[JSA[[1 ;; All,{1, 2}]] == filterSignal[[1 ;; All,{1, 2}]] == filterIdler[[1 ;; All,{1, 2}]], 
    MapThread[Append, {JSA[[1 ;; All,{1, 2}]], JSA[[1 ;; All,3]]*filterSignal[[1 ;; All,3]]*filterIdler[[1 ;; All,3]]}], 
    Print["Error: Mismatch in the PMF and PEF frequencies"]]

ComputeHeraldingEfficiency[JSA_, filterSignal_, filterIdler_] := 
    Block[{JSAfiltered, JSAfilterSignal, JSAfilterIdler}, 
    If[JSA[[1 ;; All,{1, 2}]] == filterSignal[[1 ;; All,{1, 2}]] == filterIdler[[1 ;; All,{1, 2}]], 
     JSAfilterSignal = MapThread[Append, {JSA[[1 ;; All,{1, 2}]], JSA[[1 ;; All,3]]*filterSignal[[1 ;; All,3]]}]; 
     JSAfilterIdler = MapThread[Append, {JSA[[1 ;; All,{1, 2}]], JSA[[1 ;; All,3]]*filterIdler[[1 ;; All,3]]}]; 
     JSAfiltered = ComputeFilteredJSA[JSA, filterSignal, filterIdler]; 
     Brightness[JSAfiltered]/Sqrt[Brightness[JSAfilterSignal]*Brightness[JSAfilterIdler]], 
     Print["Error: Mismatch in the PMF and PEF frequencies"]]]


OptimisePulse[frequencyMatrix_, PMF_, start_, stop_, step_, pumpFunction_, param__ ] := 
	Block[{scaling, scaledParam, tempJSA, tempPEF, purities}, 
	 scaling = Range[start, stop, step];
	 scaledParam = ReplacePart[{param},2->#*{param}[[2]]]&/@scaling;
	 tempPEF = PEFMatrix[frequencyMatrix,pumpFunction,##]&@@@ scaledParam; 
     tempJSA = JSAMatrix[#, PMF]&/@ tempPEF; 
     purities = SchmidtDecomposition[#] &/@ tempJSA; 
     Transpose[{scaling,purities[[;;,1]]}]
    ]

OptimisePulseFiltered[frequencyMatrix_, PMF_,filterSignal_, filterIdler_, start_, stop_, step_, pumpFunction_, param__ ] := 
	Block[{scaling, scaledParam, tempJSA, tempPEF, purities, tempFilteredJSA}, 
	 scaling = Range[start, stop, step];
	 scaledParam = ReplacePart[{param},2->#*{param}[[2]]]&/@scaling;
	 tempPEF = PEFMatrix[frequencyMatrix,pumpFunction,##]&@@@ scaledParam; 
     tempJSA = JSAMatrix[#, PMF]&/@ tempPEF; 
     tempFilteredJSA = ComputeFilteredJSA[#1, filterSignal,filterIdler] &/@ tempJSA; 
     purities = SchmidtDecomposition[#] &/@ tempFilteredJSA; 
     Transpose[{scaling,purities[[;;,1]]}]
    ]


SignalIdlerHOMDip[t_, decomposedJSA_] := 
    Block[{weights, signalModes, idlerModes, timeList}, 
     weights = decomposedJSA[[2]]; 
     signalModes = decomposedJSA[[3]]; 
     idlerModes = decomposedJSA[[4]];
     (1/2)*(1 - 
        Sum[weights[[i]]*weights[[j]]*
        Total[Conjugate[signalModes[[i,1 ;; All,2]]]*idlerModes[[j,1 ;; All,2]]*Exp[(-I)*signalModes[[i,1 ;; All,1]]*t]]*
        Total[signalModes[[i,1 ;; All,2]]*Conjugate[idlerModes[[j,1 ;; All,2]]]*Exp[I*idlerModes[[i,1 ;; All,1]]*t]], 
        {i, 1, Length[signalModes]}, {j, 1, Length[idlerModes]}])
        ]
SignalIdlerHOMDip[t_List, decomposedJSA_] := {#, SignalIdlerHOMDip[#, decomposedJSA]}  & /@ t

HeraldedHOMDip[t_, decomposedJSA1_, index1_, decomposedJSA2_, index2_] := 
    Block[{weights1, weights2, photon1Modes, photon2Modes,rules}, 
     rules={"s"->3,"i"->4};
     weights1 = decomposedJSA1[[2]]; 
     weights2 = decomposedJSA2[[2]]; 
     photon1Modes = decomposedJSA1[[index1/.rules]]; 
     photon2Modes = decomposedJSA2[[index2/.rules]]; 
     (1/2)*(1 - 
        Sum[weights1[[i]]^2*weights2[[j]]^2*
        Total[Conjugate[photon1Modes[[i,1 ;; All,2]]]*photon2Modes[[j,1 ;; All,2]]*Exp[(-I)*photon1Modes[[i,1 ;; All,1]]*t]]*
        Total[photon1Modes[[i,1 ;; All,2]]*Conjugate[photon2Modes[[j,1 ;; All,2]]]*Exp[I*photon2Modes[[i,1 ;; All,1]]*t]], 
        {i, 1, Length[photon1Modes]}, {j, 1, Length[photon2Modes]}])
        ]
HeraldedHOMDip[t_List, decomposedJSA1_, index1_, decomposedJSA2_, index2_] := {#, HeraldedHOMDip[#,  decomposedJSA1, index1, decomposedJSA2, index2]}  & /@ t


End[]
EndPackage[] 
