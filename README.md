# Frenkel-Kontorova+CDW Free Energy Project

## Overview

This project computes the effective free energy **F(φ,s)** of a one‑dimensional Frenkel–Kontorova (FK) chain subject to an additional charge‑density‑wave (CDW) potential. It captures the energetics of kink/anti‑kink excitations under varying CDW strength φ. Key features:

* **Minimum‑Energy Path (MEP)** via the Nudged Elastic Band (NEB) method.
* **Reaction coordinate s** defined as the cumulative arc length along the MEP.
* **Coarse‑to‑Fine φ‑scan**: accelerated annealing over φ using OpenMP parallel segments.
* **Outputs**: individual MEP CSV files, a `phi_list.csv` mapping file IDs to φ, and a combined `free_energy_surface.csv` with (φ,s,E).

---

## Directory Structure

```
project-root/
├── include/                   # Public headers
│   ├── Parameters.hpp
│   ├── InteractionPath.hpp
│   ├── SystemPhysics.hpp
│   ├── NudgedElasticBand.hpp
│   └── FreeEnergyConstruction.hpp
├── src/                       # Source files
│   ├── Parameters.cpp
│   ├── InteractionPath.cpp
│   ├── SystemPhysics.cpp
│   ├── NudgedElasticBand.cpp
│   ├── FreeEnergyConstruction.cpp
│   └── main.cpp
├── Makefile                   # Build script via GNU Make
├── README.md                  # This document
└── parameters.json            # Optional runtime parameters
```

---

## Build Instructions

Requires a C++17 compiler with OpenMP.

*Note: This project depends on two header-only libraries:*  
- **nlohmann/json** for JSON parsing  
- **Eigen** for linear algebra

### Installation

1. Clone this repo.
2. Download the JSON header and place it in `include/nlohmann/`:
   ```bash
   curl -L -o include/nlohmann/json.hpp \
     https://github.com/nlohmann/json/releases/latest/download/json.hpp
   ```
3. Download Eigen and place its `Eigen` folder in your include path:
   ```bash
   git clone https://gitlab.com/libeigen/eigen.git external/eigen
   ```
4. In your Makerfile, modify the following:
    ```make
    CXXFLAGS = -std=c++17 -O3 -fopenmp
      -I include
      -I path_to_Eigen/
      -I path_to_nlohmann_json_include/
    ```
6. Build with:
   ```bash
   make
   ```

---

This produces the executable `free_energy` in the project root.

---

## Usage

```bash
./free_energy
```

* If `parameters.json` exists, it is parsed to override defaults.
* Outputs:

  * `MEP_phi_<id>.csv` files containing atomic positions for each image at φ.
  * `phi_list.csv` mapping numeric file IDs to φ values.
  * `free_energy_surface.csv` with columns $(\varphi, s, E)$.

---

## Core Algorithm

### 1. Model Hamiltonian

The FK+CDW chain energy per configuration $X = \{x_i\}$ is:

$$
E(\varphi,X) = \sum_{i} \Bigl[ \tfrac12 K (x_{i+1}-x_i - a_0)^2 - U_0\cos(2\pi x_i / a_0) + \varphi\sin(q_{\rm cdw} x_i + \theta) \Bigr],
$$
where $\varphi$ is the CDW order parameter.

### 2. Reaction Coordinate $s$

Define $s(j) = ∑_{k=1}^j ‖X_k − X_{k−1}‖$, the cumulative Euclidean distance along an image path of M configurations.

### 3. NEB Method

* **Initialization**: linear interpolation between flat and single kink endpoints.
* **Force decomposition** into spring force along tangent and true force perpendicular to tangent.
* **Relaxation**: iterative velocity‑less updates until maximum force < tol.

### 4. Coarse‑to‑Fine $\varphi$ Scanning

1. **Coarse Sampling**: $\varphi_{\rm max} → \varphi_{\rm min}$ in steps of $\Delta\varphi_{\rm coarse} = N_{\rm threads} · \Delta\varphi$, serially.
2. **Parallel Refinement**: each coarse segment refines internal $\varphi$ points ($\Delta\varphi$) in parallel threads, using the previous image path as initial guess.
3. **Continuity Guarantee**: each segment’s initial MEP ensures physical continuity and rapid convergence.

### 5. Data Aggregation

* **phi\_list.csv**: two columns `id,phi`; `id` matches MEP file suffix.
* **free\_energy\_surface.csv**: sorted by φ descending; each row $(\varphi, s_j, E_j)$ for each image $j$.

---

## Tuning & Extensions

* **Parameters**: edit `parameters.json` or override in code.
* **Climbing‑Image NEB**: can be added for more accurate saddle‑point location.
* **Parallelization**: further via MPI or GPU for large-scale systems.

---

For questions or contributions, please open an issue or pull request on the repository.
