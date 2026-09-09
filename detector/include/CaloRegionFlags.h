#ifndef DETECTOR_CALOREGIONFLAGS_H
#define DETECTOR_CALOREGIONFLAGS_H 1

#include "DD4hep/Detector.h"

namespace k4geo {

/// Optional global switch for partial calorimeter construction (barrel-only /
/// endcap-only builds).
///
/// Looks up a global XML constant (e.g. "calo_construct_barrel" or
/// "calo_construct_endcap") and returns false only if it is defined AND equal
/// to zero. If the constant is not defined at all, this returns true, so
/// compact files that never define these constants (e.g. IDEA_o1_v03,
/// IDEA_o1_v04, IDEA_o2_v01) keep building both regions by default.
///
/// This is intentionally implemented as a *global* constant lookup rather
/// than a per-detector XML attribute (e.g. <barrel construct="..."/>)
/// because the geometry XML files driving these plugins
/// (FiberDualReadoutCalo_o1_v01.xml, SCEPCal_o1_v01.xml) are shared between
/// several IDEA options; adding a mandatory attribute there would force
/// every option that includes them to also define the same constant.
inline bool caloRegionEnabled(const dd4hep::Detector& description, const char* constantName) {
  const auto& constants = description.constants();
  auto it = constants.find(constantName);
  if (it == constants.end()) {
    return true;
  }
  return description.constantAsLong(constantName) != 0;
}

} // namespace k4geo

#endif // DETECTOR_CALOREGIONFLAGS_H
