#include "system_physics.hpp"


SystemPhysics::SystemPhysics(const Parameters& pars) : p(pars) {}

// Energy of every image (column‑wise vector)
Eigen::VectorXd SystemPhysics::energy(const InteractionPath& path, double phi) const {
    const int N = p.nAtoms;
    const int M = p.nImages;
    const double L = 2.0 * M_PI / p.q_cdw;                 // CDW period

    const Eigen::MatrixXd& X = path.X;                     // alias
    // forward neighbour under PBC
    Eigen::MatrixXd Xip = X;
    Xip.topRows(N-1) = X.bottomRows(N-1);
    Xip.row(N-1)     = X.row(0);

    // raw spring displacement Δx = x_{i+1} − x_i − a_eq
    Eigen::MatrixXd dx = Xip - X - Eigen::MatrixXd::Constant(N, M, p.a_eq);

    // minimum‑image wrap:  Δx -= L * floor( Δx/L + 0.5 )
    dx.array() -= L * ( (dx.array() / L + 0.5).floor() );

    Eigen::VectorXd E_spring = 0.5 * p.K * dx.colwise().squaredNorm();

    Eigen::MatrixXd argSub = (2.0 * M_PI / p.a_eq) * X;
    Eigen::VectorXd E_sub  = -p.U0 * argSub.array().cos().colwise().sum();

    Eigen::MatrixXd argCDW = p.q_cdw * X + Eigen::MatrixXd::Constant(N, M, p.theta);
    Eigen::VectorXd E_cdw  =  phi * argCDW.array().sin().colwise().sum();

    return E_spring + E_sub + E_cdw;
}

// True force of every image (each column vector is the true force o each image)
Eigen::MatrixXd SystemPhysics::force(const InteractionPath& path, double phi) const {
    const int N = p.nAtoms;
    const int M = p.nImages;
    const double L = 2.0 * M_PI / p.q_cdw;

    const Eigen::MatrixXd& X = path.X;

    // neighbours i‑1, i+1 under PBC
    Eigen::MatrixXd Xm = X;
    Xm.bottomRows(N-1) = X.topRows(N-1);
    Xm.row(0) = X.row(N-1);
    Eigen::MatrixXd Xp = X;
    Xp.topRows(N-1) = X.bottomRows(N-1);
    Xp.row(N-1) = X.row(0);

    Eigen::MatrixXd df = Xp - X - Eigen::MatrixXd::Constant(N, M, p.a_eq);
    Eigen::MatrixXd db = X  - Xm - Eigen::MatrixXd::Constant(N, M, p.a_eq);

    // wrap both df, db
    df.array() -= L * ((df.array()/L + 0.5).floor());
    db.array() -= L * ((db.array()/L + 0.5).floor());

    Eigen::MatrixXd grad = -p.K * (df - db);
    grad.noalias() += p.U0 * (2.0*M_PI/p.a_eq) * ((2.0*M_PI/p.a_eq)*X).array().sin().matrix();
    grad.noalias() += phi * p.q_cdw * (p.q_cdw*X + Eigen::MatrixXd::Constant(N,M,p.theta)).array().cos().matrix();
    return -grad;
}