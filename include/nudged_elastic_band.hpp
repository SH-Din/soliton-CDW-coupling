#pragma once
#include "parameters.hpp"
#include "system_physics.hpp"
#include "interaction_path.hpp"


// Class: NudgedElasticBand
// nudged elastic band (NEB) method for minimum energy path (MEP) searching
class NudgedElasticBand {
    const Parameters &p;
    const SystemPhysics &phy;
public:
    NudgedElasticBand(const Parameters &par, const SystemPhysics &physics);
    void relaxPath(InteractionPath &path, double phi);
};