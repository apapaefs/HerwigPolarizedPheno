#include "PolJetShapesHard.hh"
#include <cassert>
using namespace Rivet;
using namespace Rivet::PolJetShapesHard;
int main() {
  assert(ptWindow(19.999) == -1);
  assert(ptWindow(20) == 0);
  assert(ptWindow(29.999) == 0);
  assert(ptWindow(30) == 1);
  assert(ptWindow(44.999) == 1);
  assert(ptWindow(45) == -1);
  Vector3 normal;
  assert(!planeNormal(Vector3(0,0,1), Vector3(0,0,2), normal));
  assert(planeNormal(Vector3(1,0,0), Vector3(0,1,0), normal));
  assert((normal-Vector3(0,0,1)).mod() < 1e-12);
  double angle = 0;
  assert(signedPlaneAngle(Vector3(1,0,0), Vector3(0,1,0), Vector3(0,0,1), angle));
  assert(std::abs(angle-M_PI/2) < 1e-12);
  assert(signedPlaneAngle(Vector3(1,0,0), Vector3(0,-1,0), Vector3(0,0,1), angle));
  assert(std::abs(angle+M_PI/2) < 1e-12);
  PlaneSplit a, b, selected;
  a.z=0.3; a.kt=1; b.z=0.2; b.kt=2;
  assert(highestKtSplit({a,b}, 0.1, 1, selected));
  assert(selected.kt == 2);
  assert(!highestKtSplit({a,b}, 0.1, 2, selected));
  assert(!highestKtSplit({a,b}, 0.25, 1, selected, 0.4));
  a.kt=1.5;
  assert(highestKtSplit({a,b}, 0.25, 1, selected, 0.4));
  assert(selected.kt == 1.5);
}
