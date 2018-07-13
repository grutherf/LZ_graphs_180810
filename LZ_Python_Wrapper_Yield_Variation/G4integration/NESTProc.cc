//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The NEST program is intended for use with the Geant4 software,   *
// * which is copyright of the Copyright Holders of the Geant4        *
// * Collaboration. This additional software is copyright of the NEST *
// * development team. As such, it is subject to the terms and        *
// * conditions of both the Geant4 License, included with your copy   *
// * of Geant4 and available at http://cern.ch/geant4/license, as     *
// * well as the NEST License included with the download of NEST and  *
// * available at http://nest.physics.ucdavis.edu/                    *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutions, nor the agencies providing financial support for   *
// * this work make any representation or warranty, express or        *
// * implied, regarding this software system, or assume any liability *
// * for its use. Please read the pdf license or view it online       *
// * before download for the full disclaimer and lack of liability.   *
// *                                                                  *
// * This code implementation is based on work by Peter Gumplinger    *
// * and his fellow collaborators on Geant4 and is distributed with   *
// * the express written consent of the Geant4 collaboration. By      *
// * using, copying, modifying, or sharing the software (or any work  *
// * based on the software) you agree to acknowledge use of both NEST *
// * and Geant4 in resulting scientific publications, and you         *
// * indicate your acceptance of all the terms and conditions of the  *
// * licenses, which must always be included with this code.          *
// ********************************************************************
//
//
////////////////////////////////////////////////////////////////////////

#include "G4ParticleTypes.hh" //lets you refer to G4OpticalPhoton, etc.
#include "G4EmProcessSubType.hh" //lets you call this process Scintillation
#include "G4Version.hh" //tells you what Geant4 version you are running
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"
#include "G4ProductionCuts.hh"
#include "G4Electron.hh"
#include <cmath>
#include "NESTProc.hh"
#include "NESTStackingAction.hh"
#include "G4RandomDirection.hh"

using namespace NEST;

NESTProc<G4VUserTrackInformation> proc;

template<class T>
NESTProc<T>::NESTProc(const G4String& processName, G4ProcessType type, double efield)
      : G4VRestDiscreteProcess(processName, type), efield(efield)
{
    
        fNESTcalc = std::unique_ptr<NEST::NESTcalc>(new NEST::NESTcalc());
    
        SetProcessSubType(fScintillation);
	
        fTrackSecondariesFirst = false;
	
        if (verboseLevel>0) {
	  G4cout << GetProcessName() << " is created " << G4endl;
        }
}
template<class T>
NESTProc<T>::~NESTProc(){} //destructor needed to avoid linker error



G4Track* MakePhoton(G4ThreeVector xyz, double t) {
    // Determine polarization of new photon
    G4ParticleMomentum photonMomentum(G4RandomDirection());
    G4ThreeVector perp = photonMomentum.cross(G4RandomDirection());
    G4ThreeVector photonPolarization = perp.unit();

    G4double PhotMean = 6.97*eV; G4double PhotWidth = 0.23*eV;
    G4double sampledEnergy = G4RandGauss::shoot(PhotMean, PhotWidth);
    G4DynamicParticle* aQuantum =
            new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),
            photonMomentum);
    aQuantum->SetPolarization(photonPolarization.x(),
            photonPolarization.y(),
            photonPolarization.z());
    aQuantum->SetKineticEnergy(sampledEnergy);
    //calculate time

    return new G4Track(aQuantum,t,xyz);
    
}

template<class T>
G4VParticleChange*
NESTProc<T>::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  //ready to pop out OP and TE?
  if (NESTStackingAction::theStackingAction->isUrgentEmpty()
      && aStep.GetSecondary()->empty())
  {
    lineages_prevEvent.clear();
    photons_prevEvent.clear();

    for (auto lineage : lineages)
    {
      double etot = std::accumulate(lineage.hits.begin(), lineage.hits.end(), 0., [](double a, Hit b){return a + b.E;});
      lineage.result = fNESTcalc->FullCalculation(lineage.type, etot, lineage.density, efield, lineage.A, lineage.Z);
      auto photontimes = lineage.result.photon_times.begin();
      double ecum=0;
      double ecum_p=0;
      const double e_p = etot / lineage.result.quanta.photons;
      for (auto hit : lineage.hits)
      {
        ecum+= hit.E;
        while (ecum_p < ecum)
        {
          G4Track* onePhoton = MakePhoton(hit.xyz, *photontimes + hit.t);
          if(YieldFactor==1) aParticleChange.AddSecondary(onePhoton);
          else if (YieldFactor>0){
            if(RandomGen::rndm()->rand_uniform()<YieldFactor) aParticleChange.AddSecondary(onePhoton);
          }
          ecum_p+=e_p;
          photontimes++;
          photons_prevEvent.emplace_back(*onePhoton);
        }

      }
      assert(ecum == etot);
      lineages_prevEvent.push_back(lineage);
    }
    lineages.clear();
  }

  return G4VRestDiscreteProcess::AtRestDoIt(aTrack, aStep);

}

template<class T>
void NESTProc<T>::FillSecondaryInfo(const std::vector<G4Track*>& secondaries, NESTTrackInformation* parentInfo) const
{


    for (G4Track* sec : secondaries){
        NESTTrackInformation* infoNew = new NESTTrackInformation(*parentInfo);
        sec->SetUserInformation(infoNew);    
    }

}

template<class T>
Lineage NESTProc<T>::GetChildType(const G4Track* aTrack, const G4Track* sec) const
{
  //logic to determine what processes are kicked off by this track and also set the info

  G4String sec_creator="";
  if(sec->GetCreatorProcess()){
   sec_creator = sec->GetCreatorProcess()->GetProcessName();
  }
  if (aTrack && aTrack->GetDefinition() == G4Neutron::Definition())
  {
    return Lineage(NR);
  } else if (aTrack && aTrack->GetDefinition() == G4Gamma::Definition())
  {

    if (sec_creator.contains("compt"))
    {
      return Lineage(beta);
    } else if (sec_creator.contains("phot"))
    {
      return Lineage(gammaRay);
    }
  } else if (sec->GetDefinition() == G4Electron::Definition() && (sec_creator.contains("decay")|| !aTrack))
  {
    return Lineage(beta);
  } else if (sec->GetDefinition()->GetAtomicMass()>1 && (sec_creator.contains("decay") || !aTrack)){
    Lineage ion_lin = Lineage(ion);
    ion_lin.A=sec->GetDefinition()->GetAtomicMass();
    ion_lin.Z=sec->GetDefinition()->GetAtomicNumber();
    return ion_lin;
  }
  
  return Lineage(NoneType);
}


template<class T>
G4VParticleChange* NESTProc<T>::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
//   If a user is doing other UserTrackInfo stuff, grab it and preserve it so we don't lose the information
  G4VUserTrackInformation* oldInfo=aTrack.GetUserInformation();
  T* oldInfo_T = static_cast<T*>(oldInfo);
  NESTTrackInformation* oldInfo_N = dynamic_cast<NESTTrackInformation* >(oldInfo_T);
  
  //Type of this step.
  INTERACTION_TYPE step_type=NoneType;
  
  //hacky way to grab secondaries created in this step without them becoming const. Don't abuse this by altering the secondaries besides setting UserTrackInfo!
  const std::vector<const G4Track*>* secondaries_c = (aStep.GetSecondaryInCurrentStep());
  std::vector< G4Track*> secondaries;
  const G4TrackVector* fSecondary = aStep.GetSecondary();
  G4int nSecondary = fSecondary->size();
  for (G4int i=nSecondary - secondaries_c->size(); i < nSecondary; i++)
  {
    secondaries.push_back((*fSecondary)[i]);
  }
  
  //If the current track is already in a lineage, its secondaries inherit that lineage.
  if(oldInfo_N && oldInfo_N->parentType!=NoneType ){
    FillSecondaryInfo(secondaries,oldInfo_N);
  }
  //otherwise, we may need to start a new lineage
  else{

    for ( G4Track* sec : secondaries){
      //Each secondary has a type (including the possible NoneType)
      Lineage sec_lin = GetChildType(&aTrack, sec);
      INTERACTION_TYPE sec_type= sec_lin.type;
      //The first secondary will change the step_type. Subsequent secondaries better have the same type as the first. If they don't, something is weird
      assert(sec_type==step_type || sec_type==NoneType || step_type==NoneType);
      //if this is the first secondary to have a non-None type, we've started a new lineage
      if(step_type==NoneType && sec_type !=NoneType){
        lineages.push_back(Lineage(sec_type));
      }
      step_type=sec_type;
      //If the secondary has a non-None type, it also gets a lineage ID.
      int lineage_id = (sec_type==NoneType ? -1: lineages.size()-1);
      //If there's old (non-NEST) user info to pass on, do that. In either case, make a new NESTTrackInfo for this secondary.
      if(oldInfo_T){
        sec->SetUserInformation(new NESTTrackInformation(sec_type,lineage_id,*oldInfo_T));
      }
      else{
        sec->SetUserInformation(new NESTTrackInformation(sec_type,lineage_id));
      }
    }
    
    
    //What if the parent is a primary? Give it a lineage just as if it were one of its own secondaries
    if(aTrack.GetParentID()==0){
      Lineage sec_lin= GetChildType(0, &aTrack);
      INTERACTION_TYPE sec_type = sec_lin.type;
      assert(sec_type==step_type || sec_type==NoneType || step_type==NoneType);
      if(step_type==NoneType && sec_type !=NoneType){
          lineages.push_back(sec_lin);
        }
        step_type=sec_type;
      int hit_id = (sec_type==NoneType ? -1: lineages.size()-1);
      G4Track& aTrack_nonc = const_cast<G4Track&>( aTrack);
      if(oldInfo_T){
          aTrack_nonc.SetUserInformation(new NESTTrackInformation(sec_type,hit_id,*oldInfo_T));
        }
        else{
          aTrack_nonc.SetUserInformation(new NESTTrackInformation(sec_type,hit_id));
        }
    }
    if(step_type!=NoneType){

    }
  }
  
//  If the current track is part of a lineage...
  G4VUserTrackInformation* trackInfo=aTrack.GetUserInformation();
  T* trackInfo_T = static_cast<T*>(trackInfo);
  NESTTrackInformation* trackInfo_N = dynamic_cast<NESTTrackInformation* >(trackInfo_T);
  if(!trackInfo_N || trackInfo_N->parentType==NoneType) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  Lineage* myLineage = &lineages.at(trackInfo_N->hit_id);
  //...if the step deposited energy...
  if (aStep.GetTotalEnergyDeposit() <= 0) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

  //... in a noble element...
  const G4Material* preMaterial = aStep.GetPreStepPoint()->GetMaterial();
  const G4Material* postMaterial = aStep.GetPostStepPoint()->GetMaterial();
    G4Element *ElementA = NULL, *ElementB = NULL;
    if (preMaterial) {
        const G4ElementVector* theElementVector1 = preMaterial->GetElementVector();
        ElementA = (*theElementVector1)[0];
    }
    if (postMaterial) {
        const G4ElementVector* theElementVector2 = postMaterial->GetElementVector();
        ElementB = (*theElementVector2)[0];
    }
    G4int z1, z2;
    G4bool NobleNow = false, NobleLater = false;
    if (ElementA) z1 = (G4int) (ElementA->GetZ());
    else z1 = -1;
    if (ElementB) z2 = (G4int) (ElementB->GetZ());
    else z2 = -1;
    if (z1 == 2 || z1 == 10 || z1 == 18 || z1 == 36 || z1 == 54) {
        NobleNow = true;

    } 
    if (z2 == 2 || z2 == 10 || z2 == 18 || z2 == 36 || z2 == 54) {
        NobleLater = true;
    } 

    if ( !NobleNow ) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);



    // ...retrieve the particle's position, time, attributes at both the 
    // beginning and the end of the current step along its track...
    G4StepPoint* pPreStepPoint = aStep.GetPreStepPoint();
    G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
    G4ThreeVector x1 = pPostStepPoint->GetPosition();
    G4ThreeVector x0 = pPreStepPoint->GetPosition();
    G4double evtStrt = pPreStepPoint->GetGlobalTime();
    G4double t0 = pPreStepPoint->GetLocalTime();
    G4double t1 = pPostStepPoint->GetLocalTime();

    G4double Density = preMaterial->GetDensity() / (g / cm3); 
    if(myLineage->density==-1) myLineage->density=Density;
    double step_E = aStep.GetTotalEnergyDeposit()/keV;
    
    
    //add this hit to the appropriate lineage
    Hit stepHit(step_E,t0, x0);
    myLineage->hits.push_back(stepHit);

    //the end (exiting)
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath
// ---------------
template<class T>
G4double NESTProc<T>::GetMeanFreePath(const G4Track&,
                                          G4double ,
                                          G4ForceCondition* condition)
{
        *condition = StronglyForced;
	// what this does is enforce the G4S1Light physics process as always
	// happening, so in effect scintillation is a meta-process on top of
	// any and all other energy depositions which may occur, just like the
	// original G4Scintillation (disregard DBL_MAX, this function makes the
	// mean free path zero really, not infinite)

        return DBL_MAX; //a C-defined constant
}

// GetMeanLifeTime
// ---------------
template<class T>
G4double NESTProc<T>::GetMeanLifeTime(const G4Track&,
                                          G4ForceCondition* condition)
{
        *condition = Forced;
	// this function and this condition has the same effect as the above
        return DBL_MAX;
}
