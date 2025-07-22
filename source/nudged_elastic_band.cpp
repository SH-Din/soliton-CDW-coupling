#include "nudged_elastic_band.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

NudgedElasticBand::NudgedElasticBand(const Parameters& par, const SystemPhysics& physics)
    : p(par), phy(physics) {}

// relax the input path to MEP
void NudgedElasticBand::relaxPath(InteractionPath& path, double phi) {
    for (int iter = 0; iter < p.nebIter; ++iter) {
        // compute forces
        MatrixXd trueForce = phy.force(path, phi);
        MatrixXd newX = path.X;
        double maxF = 0;
        for (int j = 1; j < p.nImages-1; ++j) {
            // Compute unit tangent vector between neighboring images
            VectorXd tau = path.X.col(j+1) - path.X.col(j-1);
            tau.normalize();

            // Compute spring force along tangent
            VectorXd springForce = p.k_spring * (
                (path.X.col(j+1) - path.X.col(j))
                - (path.X.col(j)   - path.X.col(j-1))
            );

            double fPar = trueForce.col(j).dot(tau);
            double sPar = springForce.dot(tau);

            // Determine the NEB force
            VectorXd totalForce = springForce * tau + (trueForce.col(j) - fPar * tau);

            // Update positions and record maximum force magnitude
            newX.col(j) += p.nebStep * totalForce;
            maxF = std::max(maxF, totalForce.cwiseAbs().maxCoeff());
        }
        path.X = newX;
        if (maxF < p.nebTol) break;       // convergence
    }
    return;
}