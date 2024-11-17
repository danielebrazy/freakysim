#include "Particle.h"
#include <cmath>
#include <cstdlib>

ParticleType* Particle::ParticleTypes_[MaxNumParticleTypes_] = {nullptr};

int Particle::NumParticleTypes_ = 0;

int Particle::findParticle(const char *name)
{

    int ParticleIndex{};
    for (; ParticleIndex < NumParticleTypes_; ++ParticleIndex)
    {
        if (ParticleTypes_[ParticleIndex]->getName() == name)
        {
            break;
        }
    }
    return ParticleIndex; // aggiungere stampa a schermo?? naaah
}

void Particle::AddParticleType(const char *name, double mass, int charge, double width)
{

    if (findParticle(name) == NumParticleTypes_ && NumParticleTypes_ <= MaxNumParticleTypes_)
    {

        ResonanceType *newParticle = new ResonanceType(name, mass, charge, width);
        ParticleTypes_[NumParticleTypes_] = newParticle;

        ++NumParticleTypes_;
    }
}

void Particle::setIndex(int index)
{
    index_ = index;
}

void Particle::setIndex(const char *name)
{
    index_ = findParticle(name);
}

double Particle::getMass() const
{
    return ParticleTypes_[index_]->getMass();
}

double Particle::getCharge() const
{
    return ParticleTypes_[index_]->getCharge();
}
double Particle::getEnergy() const
{
    double ParticleEnergy = sqrt(getMass() * getMass() + getPx() * getPx() + getPy() * getPy() + getPz() * getPz());
    return ParticleEnergy;
}

double Particle::InvMass(Particle &p)
{
    double ParticleInvMass = sqrt((getEnergy() + p.getEnergy()) * (getEnergy() + p.getEnergy()) - 
                                           ((getPx() + p.getPx()) * (getPx() + p.getPx()) + 
                                            (getPy() + p.getPy()) * (getPy() + p.getPy()) + 
                                            (getPz() + p.getPz()) * (getPz() + p.getPz())));
    return ParticleInvMass;
}

void Particle::setMomentum(Momentum momentum)
{
    momentum_ = momentum;
}

int Particle::Decay2body(Particle &dau1,Particle &dau2) const {
  if(getMass() == 0.0){
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  
  double massMot = getMass();
  double massDau1 = dau1.getMass();
  double massDau2 = dau2.getMass();

  if(index_> -1){ // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;

    massMot += ParticleTypes_[index_]->getWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.setMomentum({pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta)});//fare più leggibile porca madonna 
  dau2.setMomentum({-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta)});

  double energy = sqrt(momentum_.px*momentum_.px + momentum_.py*momentum_.py + momentum_.pz*momentum_.pz + massMot*massMot);

  double bx = momentum_.px/energy;
  double by = momentum_.py/energy;
  double bz = momentum_.pz/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}
void Particle::Boost(double bx, double by, double bz)
{

  double energy = getEnergy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*momentum_.px + by*momentum_.py + bz*momentum_.pz;
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  momentum_.px += gamma2*bp*bx + gamma*bx*energy;
  momentum_.py += gamma2*bp*by + gamma*by*energy;
  momentum_.pz += gamma2*bp*bz + gamma*bz*energy;
}