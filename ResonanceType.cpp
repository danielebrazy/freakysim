#include "ResonanceType.h"

void ResonanceType::print()const{
    ParticleType::print();
    std::cout<<"width: "<<width_<<'\n';
}