// Controlled on-shell events for selection/covariance regression, not predictions.
#include "COMPASSMultiplicityFiducial.hh"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/WriterAscii.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

using namespace HepMC3;
constexpr double MP=.9382720813, MPI=.13957039, MK=.493677, MMU=.1056583755;
struct V {
  double x,y,z,e;
  V operator+(const V& b) const {return {x+b.x,y+b.y,z+b.z,e+b.e};}
  V operator-(const V& b) const {return {x-b.x,y-b.y,z-b.z,e-b.e};}
  double mass2() const {return e*e-x*x-y*y-z*z;}
};
GenParticlePtr particle(V p, int pid, int status=1) {
  return std::make_shared<GenParticle>(FourVector(p.x,p.y,p.z,p.e),pid,status);
}
V boost(V p, V parent) {
  const double mass=std::sqrt(parent.mass2()), gamma=parent.e/mass;
  const V b={parent.x/parent.e,parent.y/parent.e,parent.z/parent.e,0.};
  const double b2=b.x*b.x+b.y*b.y+b.z*b.z;
  if (b2==0.) return p;
  const double dot=b.x*p.x+b.y*p.y+b.z*p.z;
  const double a=(gamma-1.)*dot/b2+gamma*p.e;
  return {p.x+a*b.x,p.y+a*b.y,p.z+a*b.z,gamma*(p.e+dot)};
}
void twoBody(const GenVertexPtr& v, V parent, double m1, int id1, double m2, int id2) {
  const double s=parent.mass2(), mass=std::sqrt(s);
  assert(mass>m1+m2);
  const double energy=(s+m1*m1-m2*m2)/(2.*mass);
  const double p=std::sqrt(energy*energy-m1*m1);
  v->add_particle_out(particle(boost({0.,0.,p,energy},parent),id1));
  v->add_particle_out(particle(boost({0.,0.,-p,mass-energy},parent),id2));
}
void write(const std::string& path, const GenVertexPtr& v,
           const GenParticlePtr& b1, const GenParticlePtr& b2) {
  auto run=std::make_shared<GenRunInfo>();
  run->set_weight_names({"nominal"});
  GenEvent event(run,Units::GEV,Units::MM);
  event.set_event_number(1);
  event.weights()={1.};
  v->add_particle_in(b1); v->add_particle_in(b2);
  event.add_vertex(v); event.set_beam_particles(b1,b2);
  V balance={0.,0.,0.,0.};
  for (const auto& p : v->particles_in()) {
    const auto& q=p->momentum(); balance=balance+V{q.px(),q.py(),q.pz(),q.e()};
  }
  for (const auto& p : v->particles_out()) {
    const auto& q=p->momentum(); balance=balance-V{q.px(),q.py(),q.pz(),q.e()};
  }
  assert(std::abs(balance.x)+std::abs(balance.y)+std::abs(balance.z)+std::abs(balance.e)<1.e-8);
  auto xs=std::make_shared<GenCrossSection>(); xs->set_cross_section(1.,0.);
  event.set_cross_section(xs);
  WriterAscii writer(path,run); writer.write_event(event); writer.close();
}
void star(double energy, bool outside, const std::string& path) {
  auto v=std::make_shared<GenVertex>();
  V remaining={0.,0.,0.,energy};
  for (int jet=0; jet<3; ++jet) {
    const double phi=outside ? (jet==0 ? 0. : (jet==1 ? 1. : -1.)*std::acos(-.75))
                             : jet*2.*std::acos(-1.)/3.;
    const double pt=outside && jet==0 ? 30. : 20.;
    const double pz=outside && jet==0 ? pt*std::sinh(1.2) : 0.;
    const V p={.5*pt*std::cos(phi),.5*pt*std::sin(phi),.5*pz,
               std::sqrt(.25*(pt*pt+pz*pz)+MPI*MPI)};
    for (int pid : {211,-211}) {v->add_particle_out(particle(p,pid)); remaining=remaining-p;}
  }
  twoBody(v,remaining,MP,2212,MP,2212);
  const double beam=energy/2., p=std::sqrt(beam*beam-MP*MP);
  write(path,v,particle({0,0,p,beam},2212,4),particle({0,0,-p,beam},2212,4));
}
void compass(double y, int pid, double theta, bool multiple, const std::string& path) {
  constexpr double energy=160., x=.0098;
  const double nu=energy*y, scattered=energy-nu, q2=2.*MP*nu*x;
  const double beamP=std::sqrt(energy*energy-MMU*MMU);
  const double outP=std::sqrt(scattered*scattered-MMU*MMU);
  const double cosine=(energy*scattered-MMU*MMU-q2/2.)/(beamP*outP);
  const V beam={0,0,beamP,energy}, target={0,0,0,MP};
  const V lepton={outP*std::sqrt(1.-cosine*cosine),0,outP*cosine,scattered};
  const V q=beam-lepton;
  auto v=std::make_shared<GenVertex>(); v->add_particle_out(particle(lepton,-13));
  V remaining=target+q;
  if (!multiple) {
    const double mass=pid==321 ? MK : MPI, eh=.23*nu;
    const double p=std::sqrt(eh*eh-mass*mass);
    V h={p*std::sin(theta),0,p*std::cos(theta),eh};
    v->add_particle_out(particle(h,pid)); remaining=remaining-h;
  } else {
    const double qnorm=std::sqrt(q.x*q.x+q.z*q.z), nx=q.x/qnorm,nz=q.z/qnorm;
    for (int i=0; i<3; ++i) {
      const double eh=(i==2 ? .1 : .23)*nu, pt=i==0 ? .2 : .5;
      const double phi=i==0 ? .8 : 2.;
      const double pl=std::sqrt(eh*eh-MPI*MPI-pt*pt);
      V h={pl*nx+pt*std::cos(phi)*nz,pt*std::sin(phi),pl*nz-pt*std::cos(phi)*nx,eh};
      v->add_particle_out(particle(h,i==2 ? -211 : 211)); remaining=remaining-h;
    }
  }
  // A neutral on-shell baryon and pi0 complete the final state. Lambda balances K+ strangeness.
  twoBody(v,remaining,pid==321 ? 1.115683 : .9395654133,pid==321 ? 3122 : 2112,.1349768,111);
  write(path,v,particle(beam,-13,4),particle(target,2212,4));
}
void boundaries() {
  using namespace Rivet;
  const auto cell=COMPASSSIDIS::pionHadronCells().front();
  COMPASSSIDIS::DISKinematics dis;
  dis.x=.0098; dis.y=.45; dis.target=FourMomentum(MP,0,0,0);
  for (double mass : {MPI,MK}) {
    const double low=std::hypot(12.,mass)/cell.zLow;
    const double high=std::hypot(40.,mass)/cell.zHigh;
    for (double nu : {low-1.e-6,low+1.e-6,high-1.e-6,high+1.e-6}) {
      dis.q=FourMomentum(nu,0,0,nu+1.);
      const bool expected=nu>low && nu<high;
      assert(COMPASSSIDIS::multiplicityDISCell(cell,dis,mass)==expected);
      // A common boost leaves the scalar target-rest-frame nu unchanged.
      const V frame={0,0,1.,std::sqrt(2.)};
      V p=boost({0,0,0,MP},frame),q=boost({0,0,nu+1.,nu},frame);
      dis.target=FourMomentum(p.e,p.x,p.y,p.z); dis.q=FourMomentum(q.e,q.x,q.y,q.z);
      assert(COMPASSSIDIS::multiplicityDISCell(cell,dis,mass)==expected);
      dis.target=FourMomentum(MP,0,0,0);
    }
  }
}
int main(int argc,char** argv) {
  if (argc!=2) return 2;
  boundaries();
  const std::string dir=argv[1];
  for (int energy : {200,510}) for (bool outside : {false,true})
    star(energy,outside,dir+"/star-"+std::to_string(energy)+(outside ? "-outside" : "-three")+".hepmc");
  for (int pid : {211,321}) {
    const std::string suffix=pid==211 ? "-pion.hepmc" : "-kaon.hepmc";
    compass(.35,pid,.02,false,dir+"/compass-reject"+suffix);
    compass(.45,pid,.005,false,dir+"/compass-pass"+suffix);
  }
  compass(.45,211,0.,true,dir+"/compass-covariance.hepmc");
  std::cout << "nu boundary/boost checks passed; wrote 9 event fixtures\n";
}
