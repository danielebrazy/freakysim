#ifndef RESONANCETYPE_H
#define RESONANCETYPE_H
#include "ParticleType.h"

class ResonanceType : public ParticleType
{
public:
    ResonanceType(const char *name, double mass, int charge, double width) : ParticleType(name, mass, charge), width_{width} {};

    double getWidth() const { return width_; }

    void print() const;

private:
    const double width_;
};

#endif