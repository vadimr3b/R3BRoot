
#include "R3BRunConfiguration.h"

#include "TG4ComposedPhysicsList.h"
#include "TG4SpecialPhysicsList.h"

#include <LHEP_BERT.hh>

#include "R3BPhysicsList.h"

//_____________________________________________________________________________
R3BRunConfiguration::R3BRunConfiguration(const TString& userGeometry,
                                             const TString& specialProcess)
  : TG4RunConfiguration(userGeometry, "emStandard", specialProcess) 
{


}

//_____________________________________________________________________________
R3BRunConfiguration::~R3BRunConfiguration()
{
/// Destructor
}

//
// protected methods
//


//_____________________________________________________________________________
G4VUserPhysicsList*  R3BRunConfiguration::CreatePhysicsList()
{
/// Override the default physics list with user defined physics list;
/// LHEP_BERT physics list should be replaced with user own physics list

  TG4ComposedPhysicsList* builder = new TG4ComposedPhysicsList();
  
  // User physics list
  G4cout << G4endl;
  G4cout << "-I- R3BRunConfiguration Adding -R3B- SpecialPhysicsList ...  " << G4endl;
  builder->AddPhysicsList(new R3BPhysicsList());
    
//  G4cout << "Adding SpecialPhysicsList " << G4endl;
//  builder->AddPhysicsList(new TG4SpecialPhysicsList(
//                                 fSpecialProcessSelection.Data()));
  
  return builder;  
}  


