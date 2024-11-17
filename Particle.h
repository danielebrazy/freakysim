#ifndef PARTICLE_H
#define PARTICLE_H

#include "ResonanceType.h"

struct Momentum{
    double px{};
    double py{};
    double pz{};
};

class Particle{

public:

    Particle(const char* name, Momentum momentum): 
    momentum_{momentum.px, momentum.py, momentum.pz}, index_{findParticle(name)}
    {}

    Particle():momentum_{0., 0., 0.}, index_{-1}
    {};
    

    int getIndex(){return index_;}

    static void AddParticleType(const char* name, double mass, int charge, double width=0);

    void setIndex(int index);

    void setIndex(const char* name);

    static void printParticleTypes();

    void printParticleData();

    double getPx() const {return momentum_.px;}

    double getPy() const {return momentum_.py;}

    double getPz() const {return momentum_.pz;}

    double getCharge() const;

    double getMass() const;

    double getEnergy() const;

    double InvMass(Particle &p);

    void setMomentum(Momentum momentum);
//da virtuale:
    int Decay2body(Particle &dau1,Particle &dau2) const;

private:

    static const int MaxNumParticleTypes_{10};

    static ParticleType* ParticleTypes_ [MaxNumParticleTypes_];
    
    static int NumParticleTypes_;

    int index_;

    Momentum momentum_{0., 0., 0.};
    
    static int findParticle(const char* name);
//da virtuale:
    void Boost(double bx, double by, double bz);
    
};






















#endif