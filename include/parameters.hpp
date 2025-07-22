#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <nlohmann/json.hpp>


// Parameters: collects tunable knobs and *optionally* loads them from JSON.
// Any missing key in the JSON keeps its default value.
class Parameters {
public:
    // path dimensions
    int    nImages   {50};
    int    nAtoms    {100};

    // φ‑scan
    double phiMin    {0.0};
    double phiMax    {1.0};
    double dPhi      {0.01};

    // NEB
    int    nebIter   {400};
    double nebTol    {1e-4};
    double nebStep   {0.1};
    double k_spring  {1.0};

    // string method
    int    relaxIter {30};
    double relaxTol  {1e-6};
    double relaxAlpha{0.01};

    // FK + CDW physics
    double K     {1.2};
    double a_eq  {1.0};
    double U0    {0.05};
    double q_cdw {0.1};
    double theta {0.0};

    // default ctor OK
    Parameters() = default;

    // JSON‑file ctor : reads "parameters.json"‑style file.
    // Only fields present in JSON override defaults.
    explicit Parameters(const std::string& jsonFile) {          // start with defaults
        // read file → json object
        std::ifstream in(jsonFile);
        if (!in) throw std::runtime_error("Parameters::Parameters(): Cannot open JSON " + jsonFile);
        nlohmann::json jsonObj; in >> jsonObj;

        // helper macro to patch a member if key exists
#define PATCH(key, member) if (jsonObj.contains(key)) jsonObj.at(key).get_to(member)

        PATCH("nImages",    nImages);
        PATCH("nAtoms",     nAtoms);
        PATCH("phiMin",     phiMin);
        PATCH("phiMax",     phiMax);
        PATCH("dPhi",       dPhi);
        PATCH("nebIter",    nebIter);
        PATCH("nebTol",     nebTol);
        PATCH("nebStep",    nebStep);
        PATCH("k_spring",   k_spring);
        PATCH("relaxIter",  relaxIter);
        PATCH("relaxTol",   relaxTol);
        PATCH("relaxAlpha", relaxAlpha);
        PATCH("K",          K);
        PATCH("a_eq",       a_eq);
        PATCH("U0",         U0);
        PATCH("q_cdw",      q_cdw);
        PATCH("theta",      theta);

#undef PATCH
    }

    // pint(): dump all parameters to the given stream (default std::cout)
    void print(std::ostream& os = std::cout) const {
        os << "nImages      = " << nImages   << "\n"
           << "nAtoms       = " << nAtoms    << "\n"
           << "phiMin       = " << phiMin    << "\n"
           << "phiMax       = " << phiMax    << "\n"
           << "dPhi         = " << dPhi      << "\n"
           << "nebIter      = " << nebIter   << "\n"
           << "nebTol       = " << nebTol    << "\n"
           << "nebStep      = " << nebStep   << "\n"
           << "k_spring     = " << k_spring  << "\n"
           << "relaxIter    = " << relaxIter << "\n"
           << "relaxTol     = " << relaxTol  << "\n"
           << "relaxAlpha   = " << relaxAlpha<< "\n"
           << "K            = " << K         << "\n"
           << "a_eq         = " << a_eq      << "\n"
           << "U0           = " << U0        << "\n"
           << "q_cdw        = " << q_cdw     << "\n"
           << "theta        = " << theta     << std::endl;
    }
};