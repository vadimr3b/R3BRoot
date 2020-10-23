#ifndef R3BNEULANDREACTION_H
#define R3BNEULANDREACTION_H

#include "TGraph.h"
#include "TVector3.h"

#include <array>
#include <string>
#include <vector>

constexpr int AZ2PDGCode(int z, int a) { return 1000000000 + z * 10000 + a * 10; }
constexpr int PDGCode2A(int pdg) { return (pdg < 1000000000 ? 0 : ((pdg / 10) % 1000)); }
constexpr int PDGCode2Z(int pdg) { return (pdg < 1000000000 ? 0 : ((pdg / 10000) % 1000)); }
constexpr std::array<int, 2> PDGCode2AZ(int pdg) { return { PDGCode2A(pdg), PDGCode2Z(pdg) }; }

constexpr int PDG_Neutron = 2112;
constexpr int PDG_Proton = 2212;
constexpr int PDG_Deuteron = AZ2PDGCode(1, 2);
constexpr int PDG_Triton = AZ2PDGCode(1, 3);
constexpr int PDG_He3 = AZ2PDGCode(2, 3);
constexpr int PDG_He4 = AZ2PDGCode(2, 4);
constexpr int PDG_Gamma = 22;
constexpr int PDG_Electron = 11;
constexpr int PDG_Positron = -11;
constexpr int PDG_C12 = AZ2PDGCode(6, 12);
constexpr int PDG_C13 = AZ2PDGCode(6, 13);
constexpr int PDG_Al27 = AZ2PDGCode(13, 27);
constexpr int PDG_PiPlus = 211;
constexpr int PDG_PiMinus = -PDG_PiPlus;
constexpr int PDG_Pi0 = 111;

class R3BNeulandReaction
{
  public:
    static TGraph LightOutput_1H;
    static TGraph LightOutput_2H;
    static TGraph LightOutput_3He;
    static TGraph LightOutput_4He;
    static TGraph LightOutput_12C;

    struct Particle
    {
        Particle() = default;
        Particle(int pdgCode, double kineticEnergy, const TVector3& momentum)
        {
            PDGCode = pdgCode;
            KineticEnergy = kineticEnergy;
            Momentum = momentum;
        }
        int PDGCode;
        double KineticEnergy;
        TVector3 Momentum;
    };

    R3BNeulandReaction();
    virtual ~R3BNeulandReaction() {}

    void SetBarPart(int part) { fPart = part; }
    int GetBarPart() { return fPart; }

    void SetInitialNeutronEnergy(double energy) { fInitialNeutronEnergy = energy; }
    double GetInitialNeutronEnergy() { return fInitialNeutronEnergy; }

    double GetRemainingNeutronEnergy() { return fRemainingNeutronEnergy; }

    int GetReactionTyp() { return fReactionTyp; }
    int GetReactionPartner() { return fReactionPartner; }

    void SetVertex(double x, double y, double z) { fVertex = { x, y, z }; }
    void SetVertex(const TVector3& vertex) { fVertex = vertex; }
    TVector3 GetVertex() { return fVertex; };

    TVector3 GetOutgoingNeutronMomentum() { return fOutgoingNeutronMomentum; };

    void SetInteractionTime(double time) { fInteractionTime = time; }
    double GetInteractionTime() { return fInteractionTime; }

    void AddReactionProduct(int pdgCode, double kineticEnergy, const TVector3& momentum)
    {
        fParticles.emplace_back(pdgCode, kineticEnergy, momentum);
    }

    std::vector<Particle>& GetParticles() { return fParticles; }

    void Finalize();

    void Print();

    void Reset()
    {
        fParticles.clear();
        fPart = -1;
        fReactionTyp = -1;
        fRemainingNeutronEnergy = 0.;
        fOutgoingNeutronMomentum = { 0., 0., 0. };
        fMaxLightProduced = 0.;
        fMaxTrackLength = 0.;
        fInteractionTime = 1e3;
        fReactionPartner = 0;
        fVertex = { 0., 0., 0. };
    }

    static double GetMaxLight(int PDGCode, double kineticEnergy);

  private:
    std::string PDG2String(int pdgcode);
    std::string Part2String(int part);

    int fPart;
    int fReactionTyp;
    int fReactionPartner;
    double fInitialNeutronEnergy;
    double fRemainingNeutronEnergy;
    TVector3 fOutgoingNeutronMomentum;
    double fMaxLightProduced;
    double fMaxTrackLength;
    double fInteractionTime;
    TVector3 fVertex;
    std::vector<Particle> fParticles;
    ClassDef(R3BNeulandReaction, 1)
};

#endif
