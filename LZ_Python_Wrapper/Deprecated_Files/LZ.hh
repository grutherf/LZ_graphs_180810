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
double g1 = .075; //phd per S1 photon in liquid at dtCntr (not phe)
double sp = 0.38; //single phe resolution (Gaussian assumed)
double ee = 0.95; //e- extraction efficiency, overriding estimate from field
double g1_gas = 0.07; //phd per S2 photon in gas, used to get SE size
double ff = 2.0; //Fano-like fudge factor for SE width

//S1 PDE quartic poly for function of drift
//s1polA + s1polB*drift + s1polC*drift^2+... (QE included, for binomial)
double s1polA = 0.066541; //will include fit to Slava here
double s1polB =2.9440e-5;
double s1polC=-1.0780e-7;
double s1polD =3.1961e-10;
double s1polE=-2.0117e-13;

double eDriftSpeed = 1.9; // mm per usec.
double liquidBorder = 1481.; // mm

//Electric field EF as function of drift time in usec.
//The coefficients for a quintic polynomial in rising order
double efpolA = 300.; // V/cm
double efpolB = 0.;
double efpolC = 0.;
double efpolD = 0.;
double efpolE = 0.;
double efpolF = 0.;

double af = 10600.; //field in V/cm between liquid/gas border and anode

int nFold = 3; //how many PMTs have to fire for an S1 to count
double spEff = 0.0; //POD threshold in phe
double spEff_direct = 0.90; //actual efficiency, used in lieu of POD threshold

//detector parameters #1
double dtCntr = 380.; //center of detector for S1 corrections, in usec.
double dt_min = 95.0; //minimum. Top of detector fiducial
double dt_max = 750.; //maximum. Bottom of detector fiducial
double e_life = 750.; //the electron lifetime in usec.
double s2_thr = 250.; //the S2 threshold in phd. Effects NR most
double s2_rat = 0.40; //S2 bottom/total ratio, not really used anymore
double gasRho = 15e-3; //gas density in g/cm^3 based on T and P of detector
double gasGap = 0.60; //EL gap in cm, effecting both field and linear S2 term
double elGain = 0.140; //slope term for S2 photons per electron formula
double elThld = 0.474; //intercept term, setting minimum field to see S2

//detector parameters #2
int numberPMTs = 488; //in LUX, 122 minus 3 de-activated. For coincidence
double erf_mu = 1.61; //for empirical efficiency. Currently using nFold, spEff
double erf_sg = 0.61; //ditto, unused
double Xe_rho = 2.89; //for slight adjustment NEST uses liquid density. g/cm^3
double P_dphe = 0.2; //chance 1 photon makes 2 phe instead of 1 in PMT
double logMax = 2.5; //maximum log(S2/S1) admittied into analysis for limit
double logMin = 0.5; //minimum. As above, currently ignored

double vuv1_1to1_2 = 1.; //correction between different gain versions in LUX

double coinWind = 100.; //S1 concidence window in ns

//discrimination test values
double s1_min = 0.0;
double s1_max = 51.; //units of detected photons [phd]
