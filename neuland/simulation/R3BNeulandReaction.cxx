#include "R3BNeulandReaction.h"

#include "FairLogger.h"

#include <algorithm>

R3BNeulandReaction::R3BNeulandReaction()
{
    fParticles.reserve(32);
    Reset();
}

enum class ReactionTypes
{
    CarbonBreakup,
    ElasticScattering,
    ElasticProtonScattering,
    QuasiFreeProtonScattering,
    N2N,
    Unknown
};

std::string ReactionTypesString[] = {
    "Carbon Breakup", "Elastic Scattering", "Elastic Proton Scattering", "Quasi-Free Proton Scattering", "N2N",
    "Unknown"
};

void R3BNeulandReaction::Finalize()
{
    if (fParticles.size() == 0)
        return;

    fReactionTyp = static_cast<int>(ReactionTypes::Unknown);

    fRemainingNeutronEnergy = 0.;
    int A = 0;
    int Z = 0;

    auto outNPosition = fParticles.end();
    auto particleIter = fParticles.begin();
    while (true)
    {
        particleIter = std::find_if(
            particleIter, fParticles.end(), [](const Particle& particle) { return particle.PDGCode == PDG_Neutron; });
        if (particleIter == fParticles.end())
            break;

        if ((*particleIter).PDGCode == PDG_Neutron && fRemainingNeutronEnergy < (*particleIter).KineticEnergy)
        {
            outNPosition = particleIter;
            fRemainingNeutronEnergy = (*particleIter).KineticEnergy;
            fOutgoingNeutronMomentum = (*particleIter).Momentum;
        }

        ++particleIter;
    }

    if (outNPosition != fParticles.end())
        fParticles.erase(outNPosition);
    else
    {
        // Falls wir nur ein Teilchen haben, haben wir elastisch gestreut
        if (fParticles.size() == 1)
        {
            fReactionPartner = fParticles[0].PDGCode;
            fRemainingNeutronEnergy = fInitialNeutronEnergy - fParticles[0].KineticEnergy;
            if (fParticles[0].PDGCode == PDG_Proton)
            {
                fReactionTyp = static_cast<int>(ReactionTypes::ElasticProtonScattering);
                return;
            }
            else
            {
                fReactionTyp = static_cast<int>(ReactionTypes::ElasticScattering);
                return;
            }
        }
        else
        {
            // Wir haben noch immer kein Neutron.
            // Falls wir ein Pi- und ein Proton haben, wurde das Neutron in ein Proton umgewandelt

            auto protonIter = std::find_if(fParticles.begin(), fParticles.end(), [](const Particle& particle) {
                return particle.PDGCode == PDG_Proton;
            });
            auto piMinusIter = std::find_if(fParticles.begin(), fParticles.end(), [](const Particle& particle) {
                return particle.PDGCode == PDG_PiMinus;
            });

            if (protonIter != fParticles.end() && piMinusIter != fParticles.end())
            {
                auto protonPos = protonIter;
                // jetzt nehmen wir das Proton mit der höhsten Energie raus
                fRemainingNeutronEnergy = (*protonIter).KineticEnergy;
                while (++protonIter != fParticles.end())
                {
                    protonIter = std::find_if(protonIter, fParticles.end(), [](const Particle& particle) {
                        return particle.PDGCode == PDG_Proton;
                    });
                    if (protonIter == fParticles.end())
                        break;
                    if (fRemainingNeutronEnergy < (*protonIter).KineticEnergy)
                    {
                        fRemainingNeutronEnergy = (*protonIter).KineticEnergy;
                        protonPos = protonIter;
                    }
                }
                fParticles.erase(protonPos);
                ++Z;
            }
            else
            {
                // das Neutron ist in einem der Produkte.
                --A;
            }
        }
    }

    fMaxLightProduced = 0.;

    for (auto p = 0; p < fParticles.size(); ++p)
    {
        fMaxLightProduced += GetMaxLight(fParticles[p].PDGCode, fParticles[p].KineticEnergy);

        switch (fParticles[p].PDGCode)
        {
            case (PDG_Proton):
                ++A;
                ++Z;
                continue;
            case (PDG_Neutron):
                ++A;
                continue;
            case (PDG_Electron):
            case (PDG_Positron):
            case (PDG_Pi0):
            case (PDG_Gamma):
                continue;
            case (PDG_PiPlus):
                ++Z;
                continue;
            case (PDG_PiMinus):
                --Z;
                continue;
            default:
                if (fParticles[p].PDGCode > 1000000000)
                {
                    Z += (fParticles[p].PDGCode / 10000) % 1000;
                    A += (fParticles[p].PDGCode / 10) % 1000;
                }
                else
                {
                    LOG(ERROR) << "Unhandled Case: " << fParticles[p].PDGCode;
                }
        }
    }
    fReactionPartner = AZ2PDGCode(Z, A);

    if (fReactionPartner == 1000010010)
        fReactionPartner = PDG_Proton;

    if (fReactionPartner == PDG_C12)
        fReactionTyp = static_cast<int>(ReactionTypes::CarbonBreakup);
}

void R3BNeulandReaction::Print()
{
    LOG(INFO) << "Interaction in " << Part2String(fPart) << " with " << PDG2String(fReactionPartner) << ". "
              << fParticles.size() << " Particles created:";
    for (auto p = 0; p < fParticles.size(); ++p)
    {
        LOG(INFO) << "   " << PDG2String(fParticles[p].PDGCode) << " @ " << fParticles[p].KineticEnergy << " MeV";
    }
    LOG(INFO) << "Neutron Remaining Energy: " << fRemainingNeutronEnergy << " MeV.";
}

std::string R3BNeulandReaction::PDG2String(int pdgcode)
{
    switch (pdgcode)
    {
        case (PDG_Neutron):
            return "Neutron";
        case (PDG_Proton):
            return "Proton";
        case (PDG_Gamma):
            return "Gamma";
        case (PDG_C12):
            return "C12";
        case (PDG_C13):
            return "C13";
        case (PDG_Al27):
            return "Al27";
        case (PDG_He4):
            return "He4";
        case (PDG_He3):
            return "He3";
        case (PDG_Electron):
            return "Electron";
        case (PDG_Deuteron):
            return "Deuteron";
        case (PDG_Triton):
            return "Triton";
        case (PDG_Pi0):
            return "Pi0";
        case (PDG_PiMinus):
            return "Pi-";
        case (PDG_PiPlus):
            return "Pi+";
        default:
            return std::to_string(pdgcode);
    }
}

std::string R3BNeulandReaction::Part2String(int part)
{
    switch (part)
    {
        case (0):
            return "BC408";
        case (1):
            return "Aluminum";
        case (2):
            return "CH2-Wrapping";
        case (3):
            return "Air";
        default:
            return "None";
    }
}

TGraph R3BNeulandReaction::LightOutput_1H;
TGraph R3BNeulandReaction::LightOutput_2H;
TGraph R3BNeulandReaction::LightOutput_3He;
TGraph R3BNeulandReaction::LightOutput_4He;
TGraph R3BNeulandReaction::LightOutput_12C;

double R3BNeulandReaction::GetMaxLight(int PDGCode, double kineticEnergy)
{
    switch (PDGCode)
    {
        case (PDG_Proton):
            return LightOutput_1H.Eval(kineticEnergy) * kineticEnergy;
        case (PDG_Deuteron):
            return LightOutput_2H.Eval(kineticEnergy) * kineticEnergy;
        case (PDG_He3):
            return LightOutput_3He.Eval(kineticEnergy) * kineticEnergy;
        case (PDG_He4):
            return LightOutput_4He.Eval(kineticEnergy) * kineticEnergy;
        case (PDG_C12):
            return LightOutput_12C.Eval(kineticEnergy) * kineticEnergy;
        default:
            return 0.;
    }
}

ClassImp(R3BNeulandReaction)
