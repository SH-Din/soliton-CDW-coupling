#pragma once
#include <Eigen/Dense>
#include <cmath>
#include "parameters.hpp"
#include "interaction_path.hpp"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif


// SystemPhysics: Frenkel-Kontorova + CDW energy / force. 
class SystemPhysics {
    const Parameters &p;
public:
    explicit SystemPhysics(const Parameters &pars);
    Eigen::VectorXd energy(const InteractionPath &path, double phi) const;
    Eigen::MatrixXd force(const InteractionPath &path, double phi) const;
};