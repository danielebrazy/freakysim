#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H
#include <iostream>

class ParticleType{
    public:

    ParticleType(const char* name, double mass, int charge):name_{name}, mass_{mass}, charge_{charge}
    {};

    const char* getName() const {return name_;}
    double getMass() const {return mass_;}
    int getCharge() const {return charge_;}
    virtual void print() const;
    virtual double getWidth() const {return 0.;}

    private:

    const char* name_;
    const double mass_;
    const int charge_;

};

#endif