#ifndef R3BNEULANDENERGYTRACKER_H
#define R3BNEULANDENERGYTRACKER_H

#include "Rtypes.h"

#include <vector>

class R3BNeulandEnergyTracker
{
  public:
    R3BNeulandEnergyTracker();
    virtual ~R3BNeulandEnergyTracker() {}

    void SetParticleYield(int pdg, double energyloss, double light);

    const std::vector<double>& GetParticleEnergyLoss() const { return fParticleElosses; }

    const std::vector<double>& GetParticleLight() const { return fParticleLights; }

    void Reset()
    {
        for (auto& eloss : fParticleElosses)
            eloss = 0.;

        for (auto& light : fParticleLights)
            light = 0.;
    }

  private:
    std::vector<double> fParticleElosses;
    std::vector<double> fParticleLights;

    ClassDef(R3BNeulandEnergyTracker, 1)
};

#endif
