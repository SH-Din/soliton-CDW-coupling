#pragma once
#include <Eigen/Dense>
#include <string>
#include <fstream>


// Class: InteractionPath
// Stores atomic positions for all images in a single matrix
class InteractionPath {
public:
    // X(row=i, col=j) = position of atom i in image j
    Eigen::MatrixXd X;

    InteractionPath(int nAtoms, int nImages);
    static InteractionPath load(const std::string& file, int nAtoms, int nImages);
    void save(const std::string &file) const;
    static InteractionPath linearInterpolation(
        const std::string &initialFile,
        const std::string &finalFile,
        int nAtoms, int nImages);
    Eigen::VectorXd reactionCoordinate() const;
};