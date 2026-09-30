// Controlled on-shell Born-cell events for Rivet regression, not predictions.
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/WriterAscii.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace HepMC3;
constexpr double MP=.9382720813, MN=.9395654133, ME=.00051099895, MPI0=.1349768;
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
void twoBody(const GenVertexPtr& vertex, V parent, double mass, int pid) {
  const double s=parent.mass2(), total=std::sqrt(s);
  assert(total>mass+MPI0);
  const double energy=(s+mass*mass-MPI0*MPI0)/(2.*total);
  const double momentum=std::sqrt(energy*energy-mass*mass);
  vertex->add_particle_out(particle(boost({0.,0.,momentum,energy},parent),pid));
  vertex->add_particle_out(particle(boost({0.,0.,-momentum,total-energy},parent),111));
}
struct Point {int bin; double x,q2,weight;};
int main(int argc,char** argv) {
  if (argc!=4) {
    std::cerr << "Usage: generate target-pid points.txt output.hepmc\n";
    return 2;
  }
  const int pid=std::stoi(argv[1]);
  if (pid!=2212 && pid!=2112) return 2;
  const double mass=pid==2212 ? MP : MN, energy=27.6;
  std::ifstream input(argv[2]);
  std::vector<Point> points;
  Point point;
  double totalWeight=0.;
  while (input>>point.bin>>point.x>>point.q2>>point.weight) {
    points.push_back(point); totalWeight+=point.weight;
  }
  if(points.empty()) return 2;
  auto run=std::make_shared<GenRunInfo>();
  run->set_weight_names({"nominal"});
  WriterAscii writer(argv[3],run);
  for (const auto& row:points) {
    const double y=row.q2/(2.*mass*energy*row.x), scattered=energy*(1.-y);
    assert(y>0. && y<1. && scattered>ME);
    const double beamP=std::sqrt(energy*energy-ME*ME);
    const double outP=std::sqrt(scattered*scattered-ME*ME);
    const double cosine=(energy*scattered-ME*ME-row.q2/2.)/(beamP*outP);
    assert(cosine>=-1. && cosine<=1.);
    const V beam={0.,0.,beamP,energy}, target={0.,0.,0.,mass};
    const V lepton={outP*std::sqrt(1.-cosine*cosine),0.,outP*cosine,scattered};
    const V q=beam-lepton;
    assert(std::abs(-q.mass2()-row.q2)<1.e-10);
    auto vertex=std::make_shared<GenVertex>();
    vertex->add_particle_out(particle(lepton,-11));
    twoBody(vertex,target+q,mass,pid);
    auto b1=particle(beam,-11,4), b2=particle(target,pid,4);
    vertex->add_particle_in(b1); vertex->add_particle_in(b2);
    GenEvent event(run,Units::GEV,Units::MM);
    event.set_event_number(row.bin);
    event.weights()={row.weight};
    event.add_vertex(vertex); event.set_beam_particles(b1,b2);
    V balance=beam+target;
    for(const auto& p:vertex->particles_out()) {
      const auto& momentum=p->momentum();
      balance=balance-V{momentum.px(),momentum.py(),momentum.pz(),momentum.e()};
    }
    assert(std::abs(balance.x)+std::abs(balance.y)+std::abs(balance.z)+std::abs(balance.e)<1.e-8);
    auto xs=std::make_shared<GenCrossSection>(); xs->set_cross_section(totalWeight,0.);
    event.set_cross_section(xs);
    writer.write_event(event);
  }
  writer.close();
  std::cout << "Wrote " << points.size() << " on-shell e+ " << pid
            << " events, sumW=" << totalWeight << " pb normalization\n";
}
