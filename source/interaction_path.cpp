#include "interaction_path.hpp"



// Constructor: allocate storage
InteractionPath::InteractionPath(int nAtoms, int nImages)
    : X(nAtoms, nImages) {}

// Load path from CSV: each column is one image's flattened positions
InteractionPath InteractionPath::load(const std::string& file, int nAtoms, int nImages)
{
    InteractionPath path(nAtoms, nImages);
    std::ifstream in(file);
    if (!in) throw std::runtime_error("InteractionPath::load(): Cannot open file");

    std::string line;
    for (int j = 0; j < nImages; ++j) {
        if (!std::getline(in, line))
            throw std::runtime_error("InteractionPath::load(): Not enough lines for images");

        std::stringstream ss(line);
        std::string token;
        for (int i = 0; i < nAtoms; ++i) {
            if (!std::getline(ss, token, ','))           // 以逗号为分隔
                throw std::runtime_error("InteractionPath::load(): Not enough atoms in line");
            path.X(i, j) = std::stod(token);             // 转成 double
        }
    }
    return path;
}


// Save to CSV: each column on one line
void InteractionPath::save(const std::string& file) const {
    std::ofstream out(file);
    for (int j = 0; j < X.cols(); ++j) {
        for (int i = 0; i < X.rows(); ++i) {
            out << X(i, j);
            // comma between values on same line
            if (i + 1 < X.rows()) out << ',';
        }
        out << '\n';
    }
}

// create linear interpolation path from initial state to final state, loaded from files
InteractionPath InteractionPath::linearInterpolation(
    const std::string& initialFile,
    const std::string& finalFile,
    int nAtoms, int nImages
) {
    // Load single-image CSVs for the endpoints
    InteractionPath A = load(initialFile, nAtoms, 1);
    InteractionPath B = load(finalFile,   nAtoms, 1);

    // Allocate full path (nImages columns)
    InteractionPath path(nAtoms, nImages);

    // Linearly interpolate every image j (t ∈ [0,1])
    for (int j = 0; j < nImages; ++j) {
        double t = static_cast<double>(j) / (nImages - 1);
        path.X.col(j) = (1.0 - t) * A.X.col(0) + t * B.X.col(0);
    }
    return path;
}


// Compute cumulative arc length: s[j] = sum ||X[:,k] - X[:,k-1]|| for k=1..j
Eigen::VectorXd InteractionPath::reactionCoordinate() const {
    int M = X.cols();
    Eigen::VectorXd s = Eigen::VectorXd::Zero(M);
    for (int j = 1; j < M; ++j) {
        // Euclidean distance between image j and j-1
        s(j) = s(j-1) + (X.col(j) - X.col(j-1)).norm();
    }
    return s;
}