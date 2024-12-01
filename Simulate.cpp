#include "Particle.h"
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

// defining the mass of every particle type
double mPion{0.13957};
double mKaon{0.49367};
double mProton{0.93827};
double mKstar{0.89166};
// defining K* width
double wKstar{0.05};

void simulate()
{
    // adding particle types:
    Particle::AddParticleType("pi+", mPion, +1);
    Particle::AddParticleType("pi-", mPion, -1);
    Particle::AddParticleType("ka+", mKaon, +1);
    Particle::AddParticleType("ka-", mKaon, -1);
    Particle::AddParticleType("pr+", mProton, +1);
    Particle::AddParticleType("pr-", mProton, -1);
    Particle::AddParticleType("K*", mKstar, 0, wKstar);

    // setting the generation seed:
    gRandom->SetSeed();

    // defining the histograms:

    // TYPES OF PARTICLES:
    TH1F *TypesH = new TH1F("Types", "Particle Types", 7, 0.5, 7.5);
    // ANGLES:
    TH2F *AnglesH = new TH2F("Angles", "Angle Coordinates", 100, 0., 2 * M_PI, 1000, 0., M_PI);
    // MOMENTUM:
    TH1F *MomentumH = new TH1F("Momentum", "Momentum Distribution", 1000, 0., 6.);
    // TRASVERSAL MOMENTUM:
    TH1F *TransMomentumH = new TH1F("TransMomentum", "Trasversal Momentum Distribution", 1000, 0., 5.);
    // ENERGY:
    TH1F *EnergyH = new TH1F("Energy", "Energy Distribution", 1000, 0., 7);
    // TOTAL INVARIANT MASS:
    TH1F *InvMassH = new TH1F("InvMass", "Invariant Mass", 1000, 0., 9);
    // INVARIANT MASS BETWEEN PARTICLES WITH OPPOSITE CHARGE:
    TH1F *DiscInvMassH = new TH1F("DiscInvMass", "Discord Invariant Mass", 1000, 0., 9);
		DiscInvMassH->Sumw2();
    // INVARIANT MASS BETWEEN PARTICLES WITH SAME CHARGE:
    TH1F *ConcInvMassH = new TH1F("ConcInvMass", "Concord Invariant Mass", 1000, 0., 9);
		ConcInvMassH->Sumw2();
    // INVARIANT MASS BETWEEN PI+ AND KA-:
    TH1F *PpKmInvMassH = new TH1F("Pi+Ka-InvMass", "Pi+ and K- / Pi- and K+ Invariant Mass", 1000, 0., 9);
		PpKmInvMassH->Sumw2();
    // INVARIANT MASS BETWEEN PI+ AND KA+:
    TH1F *PpKpInvMassH = new TH1F("Pi+Ka+InvMass", "Pi+ and K+ / Pi- and K- Invariant Mass", 1000, 0., 9);
		PpKpInvMassH->Sumw2();
    // INVARIANT MASS BETWEEN PARTICLES DESCENDING FROM THE SAME K*:
    TH1F *DecayInvMass = new TH1F("DecayInvMass", "Decay particles invariant mass", 1000, 0., 2.);

    // defining the array which handles the generations of every event:
    Particle EventParticles[120];
    // event loop:
    for (int event{0}; event < 100000; ++event)
    {
        //std::cout << "SKIBIDI SIGMA, YOU'RE THE T O P SIGMA" <<'\n';
        // particle loop:
        for (int generation{0}; generation < 100; ++generation)
        {
            // generating random momentum:
            double phi = gRandom->Uniform(0., 2 * M_PI);
            double theta = gRandom->Uniform(0., M_PI);
            double rho = gRandom->Exp(1.);

            Momentum RndMomentum{rho * std::sin(theta) * std::cos(phi),
                                 rho * std::sin(theta) * std::sin(phi),
                                 rho * cos(theta)};

            Particle RndParticle;

            RndParticle.setMomentum(RndMomentum);

            // generating a particle of different type according to defined probabilities:
            double control = gRandom->Uniform(0., 1.);
            if (control < 0.4)
            {
                RndParticle.setIndex(0);
            }
            else if (control < 0.8)
            {
                RndParticle.setIndex(1);
            }
            else if (control < 0.85)
            {
                RndParticle.setIndex(2);
            }
            else if (control < 0.9)
            {
                RndParticle.setIndex(3);
            }
            else if (control < 0.945)
            {
                RndParticle.setIndex(4);
            }
            else if (control < 0.99)
            {
                RndParticle.setIndex(5);
            }
            else
            {
                RndParticle.setIndex(6);

                if (control < 0.995)
                {

                    Particle dau1("pi+", {0., 0., 0.});
                    Particle dau2("ka-", {0., 0., 0.});

                    RndParticle.Decay2body(dau1, dau2);

                    int i{100};
                    while (EventParticles[i].getIndex() != -1)
                    {
                        ++i;
                    }
                    EventParticles[i] = dau1;
                    EventParticles[i + 1] = dau2;

                    DecayInvMass->Fill(dau1.InvMass(dau2));
                }
                else
                {

                    Particle dau1("pi-", {0., 0., 0.});
                    Particle dau2("ka+", {0., 0., 0.});

                    RndParticle.Decay2body(dau1, dau2);
                    int i{100};
                    while (EventParticles[i].getIndex() != -1)
                    {
                        ++i;
                    }
                    EventParticles[i] = dau1;
                    EventParticles[i + 1] = dau2;

                    DecayInvMass->Fill(dau1.InvMass(dau2));
                }

                // filling the kinetic properties histograms for the generated particles:
            }
            TypesH->Fill(RndParticle.getIndex() + 1);
            AnglesH->Fill(phi, theta);
            MomentumH->Fill(std::sqrt(RndMomentum.px * RndMomentum.px + RndMomentum.py * RndMomentum.py + RndMomentum.pz * RndMomentum.pz));
            TransMomentumH->Fill(std::sqrt(RndMomentum.px * RndMomentum.px + RndMomentum.py * RndMomentum.py));
            EnergyH->Fill(RndParticle.getEnergy());

            // filling the EventParticles array:
            int i{0};
            while (EventParticles[i].getIndex() != -1)
            {
                ++i;
            }
            EventParticles[i] = RndParticle;
        }

        // reading loop for the EventParticles array, necessary to fill the invmass histos(feelin freaky:)))
        int p1{0};
        while (EventParticles[p1].getIndex() != -1)
        {
            if (EventParticles[p1].getIndex() != 6) // excluding kstar
            {
                int p2{p1 + 1};
                while (EventParticles[p2].getIndex() != -1)
                {
                    if (EventParticles[p2].getIndex() != 6) // excluding kstar
                    // calculating invariant mass:
                    {
                        double InvMass = EventParticles[p1].InvMass(EventParticles[p2]);
                        // every possible combination:
                        InvMassH->Fill(InvMass);
                        // opposite charge combinations:
                        if (EventParticles[p1].getCharge() * EventParticles[p2].getCharge() < 0)
                        {
                            DiscInvMassH->Fill(InvMass);
                            // pi+ and ka-:
                            if (EventParticles[p1].getMass() + EventParticles[p2].getMass() == mPion + mKaon)
                            {
                                PpKmInvMassH->Fill(InvMass);
                            }
                        }
                        // same charge combinations:
                        else if (EventParticles[p1].getCharge() * EventParticles[p2].getCharge() > 0)
                        {
                            ConcInvMassH->Fill(InvMass);
                            // pi+ and ka+:
                            if (EventParticles[p1].getMass() + EventParticles[p2].getMass() == mPion + mKaon)
                            {
                                PpKpInvMassH->Fill(InvMass);
                            }
                        }
                    
                    }
                    ++p2;
                }
            }
            ++p1;
        }

        // clearing the EventParticles array after every event:
        for (Particle &p : EventParticles)
        {
            p.setIndex(-1);
        }
    }

    std::cout << "Creating canvases for drawing.\n";
    TCanvas *types = new TCanvas("types", "particle types", 200, 10, 600, 400);
    TypesH->Draw();

    TCanvas *properties = new TCanvas("properties", "particle kinetic properties", 200, 10, 600, 400);
    properties->Divide(2, 2);
    properties->cd(1);

    AnglesH->Draw();
    properties->cd(2);
    MomentumH->Draw();
    properties->cd(3);
    TransMomentumH->Draw();
    properties->cd(4);
    EnergyH->Draw();

    TCanvas *invmass1 = new TCanvas("invmasstot", "total invariant mass", 200, 10, 600, 400);
    InvMassH->Draw();

    TCanvas *invmass2 = new TCanvas("invmassnontot", "invariant mass for same/opposite charge particles", 200, 10, 600, 400);
    invmass2->Divide(2, 2);
    invmass2->cd(1);
    DiscInvMassH->Draw();
    invmass2->cd(2);
    ConcInvMassH->Draw();
    invmass2->cd(3);
    PpKmInvMassH->Draw();
    invmass2->cd(4);
    PpKpInvMassH->Draw();

    TCanvas *invmass3 = new TCanvas("invmass3", "invariant mass for particles descending from the same K *", 200, 10, 600, 400);
    DecayInvMass->Draw();

    std::cout << "Done with canvases.\n";

    TFile *sigmagraphs = new TFile("Output.root", "RECREATE");

		TypesH->Write();
		AnglesH->Write();
		MomentumH->Write();
		TransMomentumH->Write();
		EnergyH->Write();
		InvMassH->Write();
		ConcInvMassH->Write();
		DiscInvMassH->Write();
		PpKpInvMassH->Write();
		PpKmInvMassH->Write();
    DecayInvMass->Write();
    sigmagraphs->Close();
}
