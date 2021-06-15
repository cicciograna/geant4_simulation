/// \file TEActionInitialization.hh
/// \brief Definition of the TEActionInitialization class

#ifndef TEActionInitialization_h
#define TEActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class TEDetectorConstruction;

/// Action initialization class.
///

class TEActionInitialization : public G4VUserActionInitialization
{
  public:
    TEActionInitialization();
    virtual ~TEActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif

    
