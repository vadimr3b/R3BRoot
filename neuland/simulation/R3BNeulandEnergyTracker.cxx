#include "R3BNeulandEnergyTracker.h"

#include "R3BNeulandReaction.h"

R3BNeulandEnergyTracker::R3BNeulandEnergyTracker()
{
    fParticleElosses.resize(8);
    fParticleLights.resize(8);
    Reset();
}

void R3BNeulandEnergyTracker::SetParticleYield(int pdg, double energyloss, double light)
{
    int index;
    switch (pdg)
    {
        case (PDG_Proton):
            index = 0;
            break;
        case (PDG_Deuteron):
            index = 1;
            break;
        case (PDG_Triton):
            index = 2;
            break;
        case (PDG_He3):
            index = 3;
            break;
        case (PDG_He4):
            index = 4;
            break;
        case (PDG_Pi0):
        case (PDG_PiPlus):
        case (PDG_PiMinus):
            index = 5;
            break;
        case (PDG_Electron):
        case (PDG_Positron):
            index = 6;
            break;
        default:
            index = 7;
    }

    fParticleElosses[index] += energyloss;
    fParticleLights[index] += light;
}