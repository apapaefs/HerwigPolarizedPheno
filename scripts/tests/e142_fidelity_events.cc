// Controlled E142 events for regression tests, not generator predictions.
#include "E142DIS.hh"
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

void electronDIS(double x, double q2, int targetId, int leptonId, const std::string& path) {
  constexpr double E=25.51, ME=0.00051099895;
  const double mass=targetId==2212 ? MP : .9395654133;
  const double Ep=E-q2/(2*mass*x), pb=std::sqrt(E*E-ME*ME), po=std::sqrt(Ep*Ep-ME*ME);
  const double costheta=(E*Ep-ME*ME-q2/2)/(pb*po);
  assert(Ep>ME && std::abs(costheta)<1.);
  V beam={0,0,pb,E}, target={0,0,0,mass};
  V scattered={po*std::sqrt(1-costheta*costheta),0,po*costheta,Ep};
  auto v=std::make_shared<GenVertex>();
  v->add_particle_out(particle(scattered,leptonId));
  twoBody(v,beam+target-scattered,mass,targetId,.1349768,111);
  write(path,v,particle(beam,leptonId,4),particle(target,targetId,4));
}
int main(int argc, char** argv) {
  assert(argc==2);
  std::string out=argv[1];
  const double means[]={.035,.05,.082,.124,.175,.246,.343,.466};
  for (int target : {2212,2112}) {
    for (size_t i=0; i<8; ++i) {
      electronDIS(means[i],E142DIS::q2Means[i],target,11,
        out+"/mean-"+std::to_string(target)+"-"+std::to_string(i)+".hepmc");
    }
  }
  electronDIS(.124,2.5*1.04,2212,11,out+"/nominal-wide.hepmc");
  electronDIS(.124,2.5*1.08,2212,11,out+"/wide-only.hepmc");
  electronDIS(.124,2.5*1.12,2212,11,out+"/outside.hepmc");
  electronDIS(.035,.995,2212,11,out+"/below-q2-floor.hepmc");
  electronDIS(.124,2.5,2212,-11,out+"/wrong-lepton.hepmc");
  // Compare the depolarization formula with the paper's geometric expression.
  for (double x : means) {
    const int bin=E142DIS::xBin(x);
    const double q2=E142DIS::q2Means[bin], E=25.51;
    for (double mass : {MP,.9395654133}) {
      const double nu=q2/(2*mass*x), Ep=E-nu;
      const double theta=2*std::asin(std::sqrt(q2/(4*E*Ep)));
      const double epsilon=1/(1+2*(1+nu*nu/q2)*std::pow(std::tan(theta/2),2));
      const double D=(E-epsilon*Ep)/(E*(1+epsilon*E142DIS::r1990(x,q2)));
      assert(std::abs(D-E142DIS::depolarization(x,nu/E,q2,mass))<1.e-12);
    }
  }
  assert(E142DIS::xBin(.0299)==-1 && E142DIS::xBin(.03)==0);
  assert(E142DIS::xBin(.04)==1 && E142DIS::xBin(.6)==-1);
  assert(!E142DIS::inWindow(.035,1.,.1));
}
