
#include "ResonanceType.h"

int main()
{

    ParticleType *Sigma = new ParticleType("sigma", 48., 76);
    ResonanceType *Beta = new ResonanceType("beta", 33, 22., 9.);

    ParticleType *ParticleArray[2];
    ParticleArray[0] = Sigma;
    ParticleArray[1] = Beta;

    for (int i = 0; i < 2; ++i)
    {
        ParticleArray[i]->print();
    }
    Beta->print();
}