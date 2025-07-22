#pragma once
#include "parameters.hpp"
#include "system_physics.hpp"
#include "nudged_elastic_band.hpp"

// Class: FreeEnergyConstruction
// coarse-to-fine φ scanning using OpenMP to produce F(phi,s)
class FreeEnergyConstruction {
    const Parameters &p;
    SystemPhysics phy;
    NudgedElasticBand neb;
public:
    explicit FreeEnergyConstruction(const Parameters &pars);
    void cal();
};