#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <vector>
#include <cfloat>
#include <iostream>
#include <fstream>

#define HEIGHT 168.7//cm
#define RADIUS  72.8
#define TrFOIL 0.//66802
double pmtRad[4] = {0.627,0.720,0.619,0.709}; //outer then inner
double gridsA[5] = { 0.20, 0.20, 0.20, 0.20, 0.20 };
double gridsZ[5] = { 1.00, 14.8, 160.4,161.7,167.7};
double gridsD[5] = { 75.0, 100., 75.0, 100., 75.0 }; //um
double gridsS[5] = { 0.50, 0.50, 0.50, 0.25, 0.50 };
double gridsT[5];
#define BORDER 160.9
//#define PMTliveFrac 0.548 //Al ring counted as dead
#define TOL 1e-9 //precision 1 nm

using namespace std;

typedef enum {
  GXe = 0,
  LXe = 1,
  Tef = 2,
  Qtz = 3,
} MATERIALS;

typedef enum {
  up = 1,
  dn =-1,
} DIRECTIONS;

double rand_gauss(); double QuantumEfficiency ( int array, double angle ) {
  
  if ( array == 0 ) return 0.25;//(0.27+rand_gauss()*0.05)*(1.-exp(-1.7/cos(angle)))/(1-exp(-1.7));
               else return 0.25;//(0.27+rand_gauss()*0.05)*(1.-exp(-1.7/cos(angle)))/(1-exp(-1.7));
  
}

double wavelength_nm ( double energy ) { //E, eV; lambda, nm
  
  return 6.626e-34 * 299792458./
    ( energy * 1e-9 * 1.602e-19 );

}

double direction[6]; void RaylScat ( bool nonRndm, double check );

double rand_uniform ( ) {

  return ( double ) rand ( ) / ( double ) RAND_MAX;

}

double rand_exp ( double mfp ) {
  
  return - mfp * log ( rand_uniform ( ) );
  
}

double rand_gauss ( ) {
  
  double u = rand_uniform(), v = rand_uniform();
  return sqrt(-2.*log(u))*cos(2.*M_PI*v);
  
}

double AngleToNormal(double origin[6]); bool FindWithinPMT(double PMTsize);
double dotP ( const double vec1[6], const double vec2[6] ) {
  
  return (vec1[3]-vec1[0])*(vec2[3]-vec2[0])+
         (vec1[4]-vec1[1])*(vec2[4]-vec2[1])+
         (vec1[5]-vec1[2])*(vec2[5]-vec2[2]);
  
}

double distanceSq ( double vect[6] ) {
  
  return (vect[3]-vect[0])*(vect[3]-vect[0])+
         (vect[4]-vect[1])*(vect[4]-vect[1])+
         (vect[5]-vect[2])*(vect[5]-vect[2]);
  
}

void truncation ( double level ) {
  
  double length2 = distanceSq ( direction );
  double scatter = ((level-direction[2])*sqrt(length2))/
    (direction[5]-direction[2]);
  direction[3] = direction[0] + scatter * (direction[3]-direction[0]) / sqrt(length2);
  direction[4] = direction[1] + scatter * (direction[4]-direction[1]) / sqrt(length2);
  direction[5] = direction[2] + scatter * (direction[5]-direction[2]) / sqrt(length2);
  
  return;
  
}

void reflection ( const double normal[6], bool diffuse ) {
  
  const double origin[6] =
    
    { direction[0],
      direction[1],
      direction[2],
      direction[3],
      direction[4],
      direction[5] }; int arrow;

  if ( (direction[5]-direction[2]) > 0. ) arrow = up;
                                     else arrow = dn;
  
  direction[0] = origin[3];
  direction[1] = origin[4];
  direction[2] = origin[5];

  if ( diffuse ) {
    
    //vector<double> pos(3);
  RE_TRY: RaylScat(false,0);
    //pos = randomDirection ( RADIUS );
    
    //direction[3] = pos[0];
    //direction[4] = pos[1];
    //direction[5] = pos[2];
    
    if ( fabs(sqrt(pow(direction[3],2.)+pow(direction[4],2.)) - RADIUS) > TOL ||
	 (fabs(direction[5]-direction[2]) < TOL &&
	  fabs(direction[4]-direction[1]) < TOL &&
	  fabs(direction[3]-direction[0]) < TOL) ) goto RE_TRY;
    
    if ( (direction[2] < (0.0000+TrFOIL+TOL) && direction[5] < direction[2]) ||
	 (direction[2] > (HEIGHT-TrFOIL-TOL) && direction[5] > direction[2]) ) goto RE_TRY;
    
    //if ( fabs(sqrt(pow(direction[0],2.)+pow(direction[1],2.)) - RADIUS) < TOL &&
    // ( (direction[2] < direction[5] && arrow == dn) || (direction[2] > direction[5] && arrow == up) )
    // && direction[2] > (0.0000+TrFOIL+TOL) && direction[2] < (HEIGHT-TOL-TrFOIL) ) { if(rand_uniform()<0.75) goto RE_TRY; }
    
    //if ( fabs(sqrt(pow(direction[0],2.)+pow(direction[1],2.)) - RADIUS) < TOL &&
    //if ( ((direction[2] > BORDER && direction[5] < BORDER) || (direction[2] < BORDER && direction[5] > BORDER)) && TIR ) goto RE_TRY;
    //&& direction[2] > (0.0000+TrFOIL+TOL) && direction[2] < (HEIGHT-TOL-TrFOIL) ) goto RE_TRY;
    
    return;
    
  }
  
  const double cosI = -dotP(normal, origin);
  
  const double transmit[3] = { (origin[3]-origin[0])+2.*cosI*(normal[3]-normal[0]),
			       (origin[4]-origin[1])+2.*cosI*(normal[4]-normal[1]),
			       (origin[5]-origin[2])+2.*cosI*(normal[5]-normal[2]) };
  
  double X = (sqrt((transmit[0]*transmit[0]+transmit[1]*transmit[1])*RADIUS*RADIUS-origin[3]*origin[3]*transmit[1]*transmit[1]+2.*origin[3]*origin[4]*transmit[0]*
		   transmit[1]-origin[4]*origin[4]*transmit[0]*transmit[0])-origin[3]*transmit[0]-origin[4]*transmit[1])/(transmit[0]*transmit[0]+transmit[1]*transmit[1]);
  
  direction[3] = origin[3] + X * transmit[0];
  direction[4] = origin[4] + X * transmit[1];
  direction[5] = origin[5] + X * transmit[2];
  
  return; // new dir
  
}

double refractiveIndex ( int material, double wavelength ) {
  
  switch ( material ) {
    
  case LXe:
    return 1.3805-9.5730/wavelength+12850./pow(wavelength,2.)-3.7518e6/pow(wavelength,3.)+6.4012e8/pow(wavelength,4.);
  case Qtz:
    return 1.5038+44.569/wavelength-17601./pow(wavelength,2.)+3.8406e6/pow(wavelength,3.)-1.9606e8/pow(wavelength,4.);
  case Tef:
    return 1.5;
  case GXe:
    return 1.000702;
    
  }
  
  return 1.0000; //if calling anything else
  
}

void refraction ( const double normal[6], int materials[2], double wavelength ) {
  
  double n1 = refractiveIndex ( materials[0], wavelength );
  double n2 = refractiveIndex ( materials[1], wavelength );
  
  const double magnitude = sqrt(distanceSq(direction));
  
  double origin[6] = { direction[0]/magnitude, direction[1]/magnitude, direction[2]/magnitude,
		       direction[3]/magnitude, direction[4]/magnitude, direction[5]/magnitude };
  
  const double n = n1 / n2;
  const double cosI = -dotP(normal, origin);
  const double sinT2 = n * n * ( 1. - cosI * cosI );
  if ( sinT2 > 1.0 ) {
    reflection ( normal, false ); return;
  }
  const double cosT = sqrt ( 1. - sinT2 );
  
  direction[0] = origin[3] * magnitude;
  direction[1] = origin[4] * magnitude;
  direction[2] = origin[5] * magnitude;
  
  const double transmit[3] = { n*(origin[3]-origin[0])+(n*cosI-cosT)*(normal[3]-normal[0]),
			       n*(origin[4]-origin[1])+(n*cosI-cosT)*(normal[4]-normal[1]),
			       n*(origin[5]-origin[2])+(n*cosI-cosT)*(normal[5]-normal[2]) };
  
  origin[3] *= magnitude;
  origin[4] *= magnitude;
  origin[5] *= magnitude;
  
  double X = (sqrt((transmit[0]*transmit[0]+transmit[1]*transmit[1])*RADIUS*RADIUS-origin[3]*origin[3]*transmit[1]*transmit[1]+2.*origin[3]*origin[4]*transmit[0]*
		   transmit[1]-origin[4]*origin[4]*transmit[0]*transmit[0])-origin[3]*transmit[0]-origin[4]*transmit[1])/(transmit[0]*transmit[0]+transmit[1]*transmit[1]);
  
  direction[3] = origin[3] + X * transmit[0];
  direction[4] = origin[4] + X * transmit[1];
  direction[5] = origin[5] + X * transmit[2];
  
  return; // the direction global variable has been updated
  
}

double reflectance ( int material );
double reflectivity ( int material, double wavelength,
		      int material2,const double normal[6] ) {
  
  double n1 = refractiveIndex ( material, wavelength );
  double n2 = refractiveIndex ( material2,wavelength );
  
  const double magnitude = sqrt(distanceSq(direction));
  
  double origin[6] = { direction[0]/magnitude, direction[1]/magnitude, direction[2]/magnitude,
                       direction[3]/magnitude, direction[4]/magnitude, direction[5]/magnitude };
  
  const double n = n1 / n2;
  const double cosI = -dotP(normal, origin);
  const double sinT2 = n * n * ( 1. - cosI * cosI );
  if ( sinT2 > 1.0 ) return 1.0; // total internal reflection
  const double cosT = sqrt ( 1. - sinT2 );
  
  const double r_s = ( n1 * cosI - n2 * cosT ) / ( n1 * cosI + n2 * cosT );
  const double r_p = ( n2 * cosI - n1 * cosT ) / ( n2 * cosI + n1 * cosT );
  
  return (r_s*r_s+r_p*r_p)/2.0; //averaging: orth and par
  
}

double reflectivity_Schlick ( int material, double wavelength,
			      int material2,const double normal[6] ) {
  
  double n1 = refractiveIndex ( material, wavelength );
  double n2 = refractiveIndex ( material2,wavelength );
  
  const double magnitude = sqrt(distanceSq(direction));
  
  double origin[6] = { direction[0]/magnitude, direction[1]/magnitude, direction[2]/magnitude,
                       direction[3]/magnitude, direction[4]/magnitude, direction[5]/magnitude };
  
  double r0 = ( n1 - n2 ) / ( n1 + n2 );
  r0 *= r0;
  double cosX = -dotP(normal, origin);
  
  if ( n1 > n2 ) {
    const double n = n1 / n2;
    const double sinT2 = n * n * ( 1. - cosX * cosX );
    if ( sinT2 > 1.0 ) return 1.0; // TIR
    cosX = sqrt ( 1. - sinT2 );
  }
  const double x = 1. - cosX;
  
  return r0 + ( 1. - r0 ) * x * x * x * x * x; //faster than pow func
  
}

double absorptionLength ( int material, double wavelength ) {
  
  switch ( material ) {
    
  case LXe:
    return 3e3; //in centimeters
  case Qtz:
    return 30.03082 + ( -0.301817 - 30.03082 ) / pow ( 1. + pow(wavelength/177.0802,10.5777), 0.0864376 ); //synthetic
  case Tef:
    return 0.;
  case GXe:
    return 500e2;
    
  }
  
  return 1e4; //default of "infinity" (100m)
  
}

double RaylScatLength ( int material, double wavelength ) {
  
  double temp = 732.79-13.494*wavelength+0.092174*pow(wavelength,2.)-0.00033424*pow(wavelength,3.)+6.6757e-07*pow(wavelength,4.);
  
  switch ( material ) {
    
  case LXe: return 30.;
    if ( temp < 20. ) return 20.0;
                 else return temp;
  default:
    return 500e2;
    
  } // the units are cm
  
}

int main ( int argc, char** argv ) {
  
  int k, arrow, volume, surf[2]; vector<double> pos(3); double alpha, kludge = 1.;
  long numPhotons = (long)1e6, numSurvivors[2] = { 0, 0 }, numAbsorbed[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double E, nm, length2, perp[6], absorb, r, phi, scatter, X[8], minX, check; bool firstTime;
  
  gridsT[0] = 1. - (2*gridsD[0]*1e-4) / gridsS[0];
  gridsT[1] = 1. - (2*gridsD[1]*1e-4) / gridsS[1];
  gridsT[2] = 1. - (2*gridsD[2]*1e-4) / gridsS[2];
  gridsT[3] = 1. - (2*gridsD[3]*1e-4) / gridsS[3];
  gridsT[4] = 1.;// - (1*gridsD[4]*1e-4) / gridsS[4];
  
  for ( long i = 0; i < numPhotons; i++ ) { firstTime = true; //cout << i << endl;
    E = 6.97 + 0.23 * rand_gauss(); nm = wavelength_nm(E); //totalD[0]=0.;totalD[1]=0.;//absorb[0]=rand_exp(absorptionLength(GXe,nm));absorb[1]=rand_exp(absorptionLength(LXe,nm));
    phi = 2.*M_PI*rand_uniform(); r = 68.8 * sqrt ( rand_uniform() );
    direction[0] = r * cos(phi);
    direction[1] = r * sin(phi);
    direction[2] = 16.3 + (146.9-16.3)*rand_uniform();
    //printf("init %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    if ( direction[2] < BORDER ) volume = LXe;
                            else volume = GXe;
  RAYLEIGH: check = rand_uniform();
    if ( firstTime ) { RaylScat(false,0.000); firstTime = false; }
                  else RaylScat(true ,check);
    //printf("ray2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
  BOUNCE:
    if ( (direction[5]-direction[2]) > 0. ) arrow = up;
                                       else arrow = dn;
    X[0] = direction[2] - 0.0000;
    X[1] = direction[2] - gridsZ[0];
    X[2] = direction[2] - gridsZ[1];
    X[3] = direction[2] - gridsZ[2];
    X[4] = direction[2] - BORDER;
    X[5] = direction[2] - gridsZ[3];
    X[6] = direction[2] - gridsZ[4];
    X[7] = direction[2] - HEIGHT;
    minX = 1e100; k = -1;
    for ( int j = 0; j < 8; j++ ) {
      if(fabs(X[j]) < TOL) X[j]=0.;
      if ( fabs(X[j]) < minX && fabs(X[j]) > 0. ) {
	if ( arrow == up && X[j] < 0. ) { minX = X[j]; k = j; }
	if ( arrow == dn && X[j] > 0. ) { minX = X[j]; k = j; }
      }
      else { ; }
    } alpha = AngleToNormal ( direction ); //cout << alpha << endl;
  DO_OVER:
    switch ( k ) {
    case 0:
      if ( direction[2] > (0.0000+TOL) && direction[5] < (0.0000+TOL) ) goto PMTb;
      else goto REFL;
    case 1:
      if ( ((direction[2] > gridsZ[0] && direction[5] < gridsZ[0]) || (direction[2] < gridsZ[0] && direction[5] > gridsZ[0])) && rand_uniform() > gridsT[0]*1.01 ) {
	goto BOT;
      }
      if ( arrow == up ) { k = 2; goto DO_OVER; } else { k = 0; goto DO_OVER; }
    case 2:
      if ( ((direction[2] > gridsZ[1] && direction[5] < gridsZ[1]) || (direction[2] < gridsZ[1] && direction[5] > gridsZ[1])) && rand_uniform() > gridsT[1]*1.00 ) {
	goto CTH;
      }
      if ( arrow == up ) { k = 3; goto DO_OVER; } else { k = 1; goto DO_OVER; }
    case 3:
      if ( ((direction[2] > gridsZ[2] && direction[5] < gridsZ[2]) || (direction[2] < gridsZ[2] && direction[5] > gridsZ[2])) && rand_uniform() > gridsT[2]*1.00 ) {
	goto GAT;
      }
      if ( arrow == up ) { k = 4; goto DO_OVER; } else { k = 2; goto DO_OVER; }
    case 4:
      if ( (direction[2] < BORDER && direction[5] > BORDER) || (direction[2] > BORDER && direction[5] < BORDER) ) { //cout << (180./M_PI)*alpha << endl;
	goto REFR;
      }
      if ( arrow == up ) { k = 5; goto DO_OVER; } else { k = 3; goto DO_OVER; }
    case 5:
      if ( ((direction[2] > gridsZ[3] && direction[5] < gridsZ[3]) || (direction[2] < gridsZ[3] && direction[5] > gridsZ[3])) && rand_uniform() > gridsT[3]*pow(cos(alpha/2),2) ) {
	goto ANE;
      }
      if ( arrow == up ) { k = 6; goto DO_OVER; } else { k = 4; goto DO_OVER; }
    case 6:
      if ( ((direction[2] > gridsZ[4] && direction[5] < gridsZ[4]) || (direction[2] < gridsZ[4] && direction[5] > gridsZ[4])) && rand_uniform() > gridsT[4]*1.00 ) {
	goto TOP;
      }
      if ( arrow == up ) { k = 7; goto DO_OVER; } else { k = 5; goto DO_OVER; }
    case 7:
      if ( direction[2] < (HEIGHT-TOL) && direction[5] > (HEIGHT-TOL) ) goto PMTt;
      else goto REFL;
    default:
      goto REFL;
    }
  PMTt:
    truncation ( HEIGHT ); length2 = distanceSq ( direction ); //totalD[0] += length2;
    //printf("top1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("top2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(GXe,nm));
    scatter = rand_exp(RaylScatLength(GXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[0]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    } PMTtt:
    perp[3] = direction[3];
    perp[4] = direction[4];
    perp[5] = 0.;
    perp[0] = direction[3];
    perp[1] = direction[4];
    perp[2] = 1.;
    if ( rand_uniform() < pmtRad[0] ) {
      if ( rand_uniform() < reflectivity(GXe,nm,Qtz,perp) )
	{ reflection ( perp, false ); goto BOUNCE; }
      if ( rand_uniform() < pmtRad[1] ) { //if(alpha>0.35)alpha=0.35;
	alpha = asin ( (refractiveIndex(GXe,nm)/refractiveIndex(Qtz,nm)) * sin ( alpha ) );
	if ( rand_uniform() < QuantumEfficiency(0,alpha) ) numSurvivors[0]++; else numAbsorbed[1]++;
      }
      else numAbsorbed[2]++;
    }
    else {
      if( rand_uniform() < 0.80 ) {
	direction[5] -= TrFOIL;
	reflection ( perp, true );
	//printf("top3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
	//printf("top4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
	goto BOUNCE;
      }
      else numAbsorbed[3]++;
    }
    continue;
  PMTb:
    truncation ( 0.0000 ); length2 = distanceSq ( direction ); //totalD[1] += length2;
    //printf("bot1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("bot2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(LXe,nm));
    scatter = rand_exp(RaylScatLength(LXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    } PMTbb:
    perp[3] = direction[3];
    perp[4] = direction[4];
    perp[5] = 1.;
    perp[0] = direction[3];
    perp[1] = direction[4];
    perp[2] = 0.;
    if ( rand_uniform() < pmtRad[2] ) {
      if ( rand_uniform() < reflectivity(LXe,nm,Qtz,perp) )
	{ reflection ( perp, false ); goto BOUNCE; }
      if ( rand_uniform() < pmtRad[3] ) { //if(alpha>0.35)alpha=0.35;
	alpha = asin ( (refractiveIndex(LXe,nm)/refractiveIndex(Qtz,nm)) * sin ( alpha ) );
	if ( rand_uniform() < QuantumEfficiency(1,alpha) ) numSurvivors[1]++; else numAbsorbed[5]++;
      }
      else numAbsorbed[6]++;
    }
    else {
      if( rand_uniform() < reflectance(LXe) ) {
	direction[5] += TrFOIL;
	reflection ( perp, true );
	//printf("bot3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
        //printf("bot4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
	goto BOUNCE;
      }
      else numAbsorbed[7]++;
    }
    continue;
  BOT:
    truncation ( gridsZ[0] ); length2 = distanceSq ( direction ); //totalD[1] += length2;
    //printf("bot1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("bot2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(LXe,nm));
    scatter = rand_exp(RaylScatLength(LXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[0] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[8]++; continue; }
  CTH:
    truncation ( gridsZ[1] ); length2 = distanceSq ( direction ); //totalD[1] += length2;
    //printf("cth1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("cth2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(LXe,nm));
    scatter = rand_exp(RaylScatLength(LXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[1] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[9]++; continue; }
  GAT:
    truncation ( gridsZ[2] ); length2 = distanceSq ( direction ); //totalD[1] += length2;
    //printf("gat1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("gat2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(LXe,nm));
    scatter = rand_exp(RaylScatLength(LXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[2] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[10]++; continue; }
  ANE:
    truncation ( gridsZ[3] ); length2 = distanceSq ( direction ); //totalD[0] += length2;
    //printf("ano1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("ano2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(GXe,nm));
    scatter = rand_exp(RaylScatLength(GXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[0]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[3] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[11]++; continue; }
  TOP:
    truncation ( gridsZ[4] ); length2 = distanceSq ( direction ); //totalD[0] += length2;
    //printf("gtp1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("gtp2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(GXe,nm));
    scatter = rand_exp(RaylScatLength(GXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[0]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[4] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[12]++; continue; }
  REFR:
    truncation ( BORDER ); length2 = distanceSq ( direction ); //totalD[volume] += length2;
    //printf("rer1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("rer2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(volume,nm));
    scatter = rand_exp(RaylScatLength(volume,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { if ( volume == GXe ) numAbsorbed[0]++; else numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( direction[2] < BORDER ) {
      perp[3] = direction[3];
      perp[4] = direction[4];
      perp[5] = 0.;
      perp[0] = direction[3];
      perp[1] = direction[4];
      perp[2] = 1.;
      surf[0] = LXe; surf[1] = GXe; refraction(perp, surf, nm);
      //printf("rer3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      //printf("rer4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
      if ( direction[5] > BORDER ) volume = GXe;
    }
    else {
      perp[3] = direction[3];
      perp[4] = direction[4];
      perp[5] = 1.;
      perp[0] = direction[3];
      perp[1] = direction[4];
      perp[2] = 0.;
      surf[0] = GXe; surf[1] = LXe; refraction(perp, surf, nm);
      //printf("rer3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      //printf("rer4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
      if ( direction[5] < BORDER ) volume = LXe;
    }
    goto BOUNCE;
  REFL:
    length2 = distanceSq ( direction ); //totalD[volume] += length2;
    //printf("ref1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("ref2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(volume,nm));
    scatter = rand_exp(RaylScatLength(volume,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { if ( volume == LXe ) numAbsorbed[4]++; else numAbsorbed[0]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ry1' %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < reflectance(volume) ) {
      r = sqrt(direction[3]*direction[3]+direction[4]*direction[4]);
      perp[3] = 0.;
      perp[4] = 0.;
      perp[5] = direction[5] / r;
      perp[0] = direction[3] / r;
      perp[1] = direction[4] / r;
      perp[2] = direction[5] / r;
      if(rand_uniform()<0.5)reflection ( perp, true ); else reflection ( perp, false );
      //printf("ref3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      //printf("ref4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
      goto BOUNCE;
    }
    else {
      if ( volume == GXe ) numAbsorbed[13]++; else numAbsorbed[14]++; continue;
    }
  } //end of photons loop
  
  cout << "gas\t\t\t" << numAbsorbed[0]  << endl;
  cout << "photo-cathode (T)\t" << numAbsorbed[1]  << endl;
  cout << "window (T)\t\t" << numAbsorbed[2]  << endl;
  cout << "trifoil (T)\t\t" << numAbsorbed[3]  << endl;
  cout << "liquid\t\t\t" << numAbsorbed[4]  << endl;
  cout << "photo-cathode (B)\t" << numAbsorbed[5]  << endl;
  cout << "window (B)\t\t" << numAbsorbed[6]  << endl;
  cout << "trifoil (B)\t\t" << numAbsorbed[7]  << endl;
  cout << "bot grid\t\t" << numAbsorbed[8]  << endl;
  cout << "cathode\t\t\t" << numAbsorbed[9]  << endl;
  cout << "gate\t\t\t" << numAbsorbed[10] << endl;
  cout << "anode\t\t\t" << numAbsorbed[11] << endl;
  cout << "top grid\t\t" << numAbsorbed[12] << endl;
  cout << "PTFE wall in gas\t" << numAbsorbed[13] << endl;
  cout << "PTFE wall in liquid\t" << numAbsorbed[14] << endl;
  
  cout << "#sur" << "\t" << "#abs" << "\t" << "#tot" << "\t" << "g1" << "\t\t" << "[%]" << "\t" << "TBA" << endl;
  double total = 1. * double(numSurvivors[0]) + 1. * double(numSurvivors[1]);
  long noAbs = numAbsorbed[0] + numAbsorbed[1] + numAbsorbed[2] + numAbsorbed[3] + numAbsorbed[4] + numAbsorbed[5] + numAbsorbed[6] +
    numAbsorbed[7] + numAbsorbed[8] + numAbsorbed[9] + numAbsorbed[10]+ numAbsorbed[11]+ numAbsorbed[12]+ numAbsorbed[13] +numAbsorbed[14];
  cout << total << "\t" << noAbs << "\t" << numSurvivors[0]+numSurvivors[1]+noAbs << "\t" << total
    /double(numPhotons) << "\t" << 100.*total
    /double(numPhotons) << "\t" <<
    ( 1. * double(numSurvivors[0]) - 1. * double(numSurvivors[1]) ) / total << endl;
  
  return 1;
  
}

double reflectance ( int material ) {
  
  switch ( material ) {
    
  case LXe:
    return 0.95;
  case GXe:
    return 0.20;
    
  }
  
  return 0.; //if calling anything else make it non-reflective
  
}

double AngleToNormal ( double origin[6] ) {
  
  double normal[6];
  
  normal[0] = origin[3]; normal[1] = origin[4]; normal[2] = 0.;
  normal[3] = origin[3]; normal[4] = origin[4]; normal[5] = 1.;
  
  double numer = fabs(dotP(normal,origin));
  double denom = sqrt(distanceSq(normal)*distanceSq(origin));
  
  return acos ( numer / denom ); //in radians!

}

void RaylScat ( bool nonRndm, double check ) {
  
  double CosTheta=rand_uniform(); double portion=.5;
  if ( nonRndm && check < 0.5 ) {
    double dist = rand_uniform() + 1.;
    CosTheta = sqrt ( dist - 1. ); portion = 0.499;
    //CosTheta = 1.;
    //CosTheta = 2. * rand_uniform() - 1.;
    //fcostheta = ( 1. + CosTheta*CosTheta)/2.;
  }
  double SinTheta = sqrt(1.-CosTheta*CosTheta);
  // consider for the angle 90-180 degrees
  if (rand_uniform() < portion ) CosTheta = -CosTheta; //if(nonRndm) cout << CosTheta << endl;

  // simulate the phi angle
  double rand = 2.*M_PI*rand_uniform();
  double SinPhi = sin(rand);
  double CosPhi = cos(rand);
  
  // start constructing the new momentum direction
  double unit_x = SinTheta * CosPhi;
  double unit_y = SinTheta * SinPhi;
  double unit_z = CosTheta;
  
  double alpha = (sqrt(-(direction[1]*direction[1]-RADIUS*RADIUS)*unit_x*unit_x+2.*direction[0]*direction[1]*unit_x*unit_y-(direction[0]*direction[0]-RADIUS*RADIUS)*unit_y*
		       unit_y)-direction[0]*unit_x-direction[1]*unit_y)/(unit_x*unit_x+unit_y*unit_y); //extend the unit vector until ray of light hits the wall at RADIUS
  
  direction[3] = alpha * unit_x + direction[0];
  direction[4] = alpha * unit_y + direction[1];
  direction[5] = alpha * unit_z + direction[2];
  
  return; //void because dir is global var
  
}

bool FindWithinPMT ( double PMTsize ) {
  
  double PMTy[5] = { 0., 5.19615, 10.3923, 15.5885, 20.7846 };
  double PMTx[9] = { 0., 3., 6., 9., 12., 15., 18., 21., 24 };
  
  int i, minX, minY; double minV = DBL_MAX; double sign[2];
  for ( i = 0; i < 5; i++ ) {
    if ( fabs ( direction[4] - PMTy[i] ) < minV ) { minV = fabs ( direction[4] - PMTy[i] ); minY = i; sign[1] = 1.; }
  }
  for ( i = 0; i < 5; i++ ) {
    if ( fabs ( direction[4] + PMTy[i] ) < minV ) { minV = fabs ( direction[4] + PMTy[i] ); minY = i; sign[1] =-1.; }
  }
  if ( minV > PMTsize ) return false;
  
  if ( minY == 0 ) { PMTx[1] = DBL_MAX; PMTx[3] = DBL_MAX; PMTx[5] = DBL_MAX; PMTx[7] = DBL_MAX; }
  else if ( minY == 1 ) { PMTx[0] = DBL_MAX; PMTx[2] = DBL_MAX; PMTx[4] = DBL_MAX; PMTx[6] = DBL_MAX; PMTx[8] = DBL_MAX; }
  else if ( minY == 2 ) { PMTx[1] = DBL_MAX; PMTx[3] = DBL_MAX; PMTx[5] = DBL_MAX; PMTx[7] = DBL_MAX; PMTx[8] = DBL_MAX; }
  else if ( minY == 3 ) { PMTx[0] = DBL_MAX; PMTx[2] = DBL_MAX; PMTx[4] = DBL_MAX; PMTx[6] = DBL_MAX; PMTx[7] = DBL_MAX; PMTx[8] = DBL_MAX; }
  else { PMTx[1] = DBL_MAX; PMTx[3] = DBL_MAX; PMTx[5] = DBL_MAX; PMTx[6] = DBL_MAX; PMTx[7] = DBL_MAX; PMTx[8] = DBL_MAX; }
  
  minV = DBL_MAX;
  for ( i = 0; i < 9; i++ ) {
    if ( fabs ( direction[3] - PMTx[i] ) < minV ) { minV = fabs ( direction[3] - PMTx[i] ); minX = i; sign[0] = 1.; }
  }
  for ( i = 0; i < 9; i++ ) {
    if ( fabs ( direction[3] + PMTx[i] ) < minV ) { minV = fabs ( direction[3] + PMTx[i] ); minX = i; sign[0] =-1.; }
  }
  if ( minV > PMTsize ) return false;
  
  if ( ( (direction[3]-sign[0]*PMTx[minX])*(direction[3]-sign[0]*PMTx[minX]) +
	 (direction[4]-sign[1]*PMTy[minY])*(direction[4]-sign[1]*PMTy[minY]) ) < PMTsize*PMTsize )
    return true;
  else
    return false;
  
}
