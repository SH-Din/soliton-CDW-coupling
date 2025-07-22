#include "free_energy_construction.hpp"


FreeEnergyConstruction::FreeEnergyConstruction(const Parameters &pars) : p(pars), phy(pars), neb(pars, phy) {}

// scan free energy F(phi,s), write free_energy_surface.csv
void FreeEnergyConstruction::cal() {
    // Coarse phi list
    std::vector<double> coarsePhiList;
    for (double phi = p.phiMax; phi >= p.phiMin - 1e-12; phi -= p.dPhi * omp_get_max_threads()) {
        coarsePhiList.push_back(phi);
    }
    int nCoarsePhi = coarsePhiList.size();

    // Write phi list for all saved files
    std::vector<std::pair<int, double>> phiLabels;
    int phiCnt = 0;

    std::vector<int> coarseIds(nCoarsePhi);
    InteractionPath path = InteractionPath::linearInterpolation("initial.csv", "final.csv", p.nAtoms, p.nImages);
    for (int i = 0; i < nCoarsePhi; ++i) {
        double phi = coarsePhiList[i];
        neb.relaxPath(path, phi);
        std::string fileName = "MEP_phi_" + std::to_string(phiCnt) + ".csv";
        path.save(fileName);
        phiLabels.emplace_back(phiCnt, phi);
        phiCnt ++;
    }

    // Parallel fine segments
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nCoarsePhi - 1; ++i) {
        InteractionPath segPath = InteractionPath::load("MEP_phi_" + std::to_string(coarsePhiList[i]) + ".csv", p.nAtoms, p.nImages);
        for (double phi = coarsePhiList[i] - p.dPhi; phi > coarsePhiList[i + 1] + 1e-12; phi -= p.dPhi) {
            neb.relaxPath(segPath, phi);
            int fileId;
            #pragma omp critical
            {
                fileId = phiCnt;
                phiLabels.emplace_back(fileId, phi);
                phiCnt++;
            }
            std::string fileName = "MEP_phi_" + std::to_string(fileId) + ".csv";
            segPath.save(fileName);
            
        }
    }

    // Write phi index list
    std::ofstream phiListOut("phi_list.csv");
    phiListOut << "id,phi\n";
    std::sort(phiLabels.begin(), phiLabels.end(), 
        [](auto &a, auto &b) {return a.second > b.second;});
    for (const auto &[id, phi] : phiLabels) {
        phiListOut << id << "," << std::setprecision(10) << phi << "\n";
    }

    // Write free-energy surface F(phi, s)
    std::ofstream surf("free_energy_surface.csv");
    surf << "phi,s,energy\n";
    for (const auto &[id, phi] : phiLabels) {
        std::string fn = "MEP_phi_" + std::to_string(id) + ".csv";
        InteractionPath ip = InteractionPath::load(fn, p.nAtoms, p.nImages);
        Eigen::VectorXd s = ip.reactionCoordinate();
        Eigen::VectorXd E = phy.energy(ip, phi);
        for (int j = 0; j < p.nImages; ++j) {
            surf << std::fixed << std::setprecision(8)
                    << phi << "," << s(j) << "," << E(j) << '\n';
        }
    }
}