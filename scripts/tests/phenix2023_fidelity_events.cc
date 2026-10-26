// Controlled selection/normalization fixtures, not Herwig predictions.
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/WriterAscii.h"
#include <cassert>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
using namespace HepMC3;
struct Item { double pt,eta,phi,mass; int pid; };
GenParticlePtr make(const Item& i,int status=1) {
  double px=i.pt*std::cos(i.phi),py=i.pt*std::sin(i.phi),pz=i.pt*std::sinh(i.eta);
  return std::make_shared<GenParticle>(FourVector(px,py,pz,std::sqrt(px*px+py*py+pz*pz+i.mass*i.mass)),i.pid,status);
}
GenEvent makeEvent(const std::vector<Item>& items,int ancestor=0,
                   double weight=1.,int number=1,double crossSection=1.) {
  auto run=std::make_shared<GenRunInfo>();run->set_weight_names({"nominal"});
  GenEvent ev(run,Units::GEV,Units::MM);ev.set_event_number(number);ev.weights()={weight};
  const double pb=std::sqrt(255.*255.-.9382720813*.9382720813);
  auto a=std::make_shared<GenParticle>(FourVector(0,0,pb,255),2212,4);
  auto b=std::make_shared<GenParticle>(FourVector(0,0,-pb,255),2212,4);
  auto hard=std::make_shared<GenVertex>();hard->add_particle_in(a);hard->add_particle_in(b);
  FourVector total;
  for (size_t j=0;j<items.size();++j) {
    auto p=make(items[j]);
    if (j==0 && ancestor) {
      // A soft partner gives a physical pi0 mass for the hadron-decay veto.
      const double eps=.1349768*.1349768/(4*p->momentum().e());
      auto partner=std::make_shared<GenParticle>(FourVector(-eps,0,0,eps),22,1);
      auto parent=std::make_shared<GenParticle>(p->momentum()+partner->momentum(),ancestor,2);
      auto decay=std::make_shared<GenVertex>();
      hard->add_particle_out(parent);decay->add_particle_in(parent);
      decay->add_particle_out(p);decay->add_particle_out(partner);ev.add_vertex(decay);
      total+=partner->momentum();
    } else hard->add_particle_out(p);
    total+=p->momentum();
  }
  // Balance the event with two invisible massless particles, away from any
  // isolation contribution. This keeps the selection fixtures self-contained.
  const double er=510.-total.e(),px=-total.px(),py=-total.py(),pz=-total.pz();
  const double mass=std::sqrt(er*er-px*px-py*py-pz*pz),gamma=er/mass;
  const double bx=px/er,by=py/er,bz=pz/er,b2=bx*bx+by*by+bz*bz;
  for (int sign:{-1,1}) {
    const double e=mass/2,z=sign*e,dot=bz*z;
    const double factor=(b2>0?(gamma-1)*dot/b2:0)+gamma*e;
    hard->add_particle_out(std::make_shared<GenParticle>(
      FourVector(factor*bx,factor*by,z+factor*bz,gamma*(e+dot)),sign*12,1));
  }
  ev.add_vertex(hard);ev.set_beam_particles(a,b);
  auto xs=std::make_shared<GenCrossSection>();xs->set_cross_section(crossSection,0.);ev.set_cross_section(xs);
  return ev;
}
void fixture(const std::string& file,const std::vector<Item>& items,int ancestor=0) {
  auto ev=makeEvent(items,ancestor);
  WriterAscii out(file,ev.run_info());out.write_event(ev);out.close();
}
int main(int argc,char** argv) {
  assert(argc==2);std::string dir=argv[1];
  const std::vector<double> edges={6,6.5,7,7.5,8,8.5,9,9.5,10,12,14,16,18,20,22,24,26,28,30};
  for(size_t i=0;i+1<edges.size();++i)
    fixture(dir+"/bin-"+std::to_string(i)+".hepmc",{{(edges[i]+edges[i+1])/2,0,0,0,22}});
  auto f=[&](std::string name,std::vector<Item> extra,double pt=6.2,double eta=0) {
    extra.insert(extra.begin(),{pt,eta,0,0,22});fixture(dir+"/"+name+".hepmc",extra);
  };
  f("below-pt",{},5.99);f("above-pt",{},30.01);f("outside-eta",{},6.2,.251);
  f("two-same-bin",{{7.4,0,3.141592653589793,0,22}},7.1);
  f("two-photons",{{8.,0,3.141592653589793,0,22}},7.);
  f("cone-fail",{{.8,0,.1,.13957039,211}});
  f("outside-cone",{{.8,0,.501,.13957039,211}});
  f("neutrino",{{4.,0,.1,0,14}});
  f("charged-below",{{.195,0,.1,.13957039,211},{.195,0,.1,.13957039,-211},{.4,0,.2,0,22}});
  // |p|=.205 but pT<.2: must be included with the published momentum cut.
  f("charged-momentum",{{.205/std::cosh(.4),.4,.1,.13957039,211},{.4,0,.2,0,22}});
  f("neutral-below",{{.299,0,.1,0,22},{.299,0,.1,0,22},{.299,0,.1,0,22}});
  f("neutral-above",{{.301,0,.1,0,22},{.301,0,.1,0,22},{.301,0,.1,0,22}});
  f("energy-not-et",{{.60,.4,.1,0,22}}); // E=.6487>.62 while ET=.60<.62
  // Zero accepted photons must still enter the normalization denominator.
  fixture(dir+"/zero-photons.hepmc",{});
  // Multiple weighted events, unequal photon multiplicities and a non-unit
  // generator cross section. These expose both denominator and sumW2 errors.
  WriterAscii weighted(dir+"/weighted-multiplicity.hepmc");
  weighted.write_event(makeEvent({},0,1.,1,12.));
  weighted.write_event(makeEvent({{7.1,0,0,0,22},{7.4,0,3.141592653589793,0,22}},0,2.,2,12.));
  weighted.write_event(makeEvent({{7.2,0,0,0,22}},0,3.,3,12.));
  weighted.close();
  fixture(dir+"/hadron-decay.hepmc",{{8.,0,0,0,22}},111);
  fixture(dir+"/quark-radiation.hepmc",{{8.,0,0,0,22}},2);
}
