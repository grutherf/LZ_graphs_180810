#include <iostream>
#include <cmath>
#include <fstream>

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFeldmanCousins.h"
#include "TRandom3.h"
#include "TMath.h"

#include "string.h"
#define AVO 6.022e23

//S1 (first 2) and S2 (last 3) parameters
double g1 = 0.117; //phd per S1 photon in liquid at dtCntr (not phe)
double sp = 0.37; //single phe resolution (Gaussian assumed)
double ee = 0.489; //e- extraction efficiency, overriding estimate from field
double g1_gas = 0.0946; //phd per S2 photon in gas, used to get SE size
double ff = 1.4; //Fano-like fudge factor for SE width

//S1 PDE quartic poly for function of drift
//s1polA + s1polB*drift + s1polC*drift^2+... (QE included, for binomial)
double s1polA = 0.84637;
double s1polB = 0.0016222;
double s1polC =-9.9224e-6;
double s1polD = 4.7508e-8;
double s1polE =-7.1692e-11;

double eDriftSpeed = 1.51; // mm per usec.
double liquidBorder = 544.2198; // mm

//Electric field EF as function of drift time in usec.
//The coefficients for a quintic polynomial in rising order
double efpolA=214.5734 ; // V/cm
double efpolB =-0.35537;
double efpolC = 0.00113;
double efpolD =-2.1766e-6;
double efpolE = 1.2368e-9;
double efpolF = 0.;

double af = 5650.; //field in V/cm between liquid/gas border and anode

int nFold = 2; //how many PMTs have to fire for an S1 to count
double spEff = 0.3; //POD threshold in phe
double spEff_direct = 1.0; //actual efficiency, used in lieu of POD threshold

//detector parameters #1
double dtCntr = 160.; //center of detector for S1 corrections, in usec.
double dt_min = 38.0; //minimum. Top of detector fiducial
double dt_max = 305.; //maximum. Bottom of detector fiducial
double e_life = 800.; //the electron lifetime in usec.
double s2_thr = 150.; //the S2 threshold in phd. Effects NR most
double s2_rat = 0.449; //S2 bottom/total ratio, not really used anymore
double gasRho = 15e-3; //gas density in g/cm^3 based on T and P of detector
double gasGap = 0.56; //EL gap in cm, effecting both field and linear S2 term
double elGain = 0.140; //slope term for S2 photons per electron formula
double elThld = 0.474; //intercept term, setting minimum field to see S2

//detector parameters #2
int numberPMTs = 119; //in LUX, 122 minus 3 de-activated. For coincidence
double erf_mu = 1.61; //for empirical efficiency. Currently using nFold, spEff
double erf_sg = 0.61; //ditto, unused
double Xe_rho = 2.888; //for slight adjustment NEST uses liquid density. g/cm^3
double P_dphe = 0.173; //chance 1 photon makes 2 phe instead of 1 in PMT
double logMax = 2.5; //maximum log(S2/S1) admittied into analysis for limit
double logMin = 0.5; //minimum. As above, currently ignored

double vuv1_1to1_2 = .915; //correction between different gain versions in LUX

double coinWind = 100.; //S1 concidence window in ns

//discrimination test values
double s1_min = 1.0;
double s1_max = 51.; //units of detected photons [phd]
