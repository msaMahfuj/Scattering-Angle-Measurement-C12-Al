# Carbon-Ion Scattering in Aluminum Using CR-39 Nuclear Track Detectors

This repository contains the experimental data, analysis codes, GEANT4 simulation files, figures, and results for the measurement of scattering angles of carbon-ion projectiles transmitted through an aluminum target using CR-39 nuclear track detectors.

The analysis was performed for six projectile energies:

- 21.0 MeV/u
- 28.5 MeV/u
- 83.0 MeV/u
- 102.5 MeV/u
- 107.0 MeV/u
- 110.0 MeV/u

---

## Experimental Setup

For each beam energy, two CR-39 detector sheets were placed on opposite sides of an aluminum target.

- **Sheet A:** Upstream detector
- **Sheet B:** Downstream detector

Each CR-39 sheet was scanned on both its front and back surfaces. Therefore, four datasets were analyzed for each energy:

- A-front
- A-back
- B-front
- B-back

The datasets contain the centroid coordinates of etched ion tracks extracted using PitFit. The primary coordinates used in the analysis are:

- `GravX`
- `GravY`

Coordinates are expressed in micrometers (µm).

The detector geometry is defined as:

- z_A,f = 0
- z_A,b = t_A
- z_B,f = t_A + t_Al
- z_B,b = t_A + t_Al + t_B

where:

- t_A = post-etch thickness of Sheet A
- t_B = post-etch thickness of Sheet B
- t_Al = thickness of the aluminum target

---

## Analysis Workflow

The same analysis procedure was applied independently to all six beam energies.

### 1. Data Loading

For each energy, four CSV files corresponding to the four detector surfaces are loaded:

```text
A-front
A-back
B-front
B-back
```

Entries with missing values in `GravX` or `GravY` are removed before further analysis.

---

### 2. Front-Back Alignment

The front and back surfaces of a CR-39 detector may have different coordinate origins due to scanning and measurement offsets.

For each detector sheet, the coordinate differences are calculated as:

**Δx = x_back − x_front**

**Δy = y_back − y_front**

A two-dimensional Kernel Density Estimation (KDE) is applied to the **(Δx, Δy)** distribution.

The peak of the KDE distribution gives the relative translational offset between the front and back coordinate systems.

The back-surface coordinates are then corrected using:

**x′_back = x_back − Δx_peak**

**y′_back = y_back − Δy_peak**

This alignment procedure is performed independently for both Sheet A and Sheet B.

---

### 3. Front-Back Track Matching

After alignment, tracks on the front and back surfaces of each CR-39 sheet are matched using a KD-tree nearest-neighbour search.

The matching radius is determined from the spatial spread of the alignment distribution:

**r = 4σ**

Only unique one-to-one matches are retained. Ambiguous matches and tracks without a valid counterpart are excluded.

This produces reconstructed track pairs for:

- **Sheet A:** Incident trajectories before the aluminum target.
- **Sheet B:** Outgoing trajectories after the aluminum target.

---

### 4. Three-Dimensional Trajectory Reconstruction

For every matched front-back track pair, a three-dimensional trajectory vector is reconstructed.

For Sheet A:

**v⃗_A = (x_A,b − x_A,f, y_A,b − y_A,f, dz_A)**

For Sheet B:

**v⃗_B = (x_B,b − x_B,f, y_B,b − y_B,f, dz_B)**

where:

- dz_A = post-etch thickness of Sheet A
- dz_B = post-etch thickness of Sheet B

The corresponding unit vectors are calculated as:

**v̂ = v⃗ / |v⃗|**

The trajectory reconstructed in Sheet A represents the incident direction of the carbon ion before passing through the aluminum target, while the trajectory reconstructed in Sheet B represents the outgoing direction after transmission through the target.

---

### 5. Cross-Matching of Incident and Outgoing Tracks

The trajectory reconstructed in Sheet A is projected to the plane of the front surface of Sheet B using the reconstructed direction and the known separation between the detector sheets.

The projected position is then compared with the reconstructed tracks in Sheet B using a KD-tree nearest-neighbour search.

Only successfully cross-matched trajectories are retained for the final scattering-angle calculation.

---

### 6. Scattering-Angle Calculation

For each cross-matched event, the scattering angle is calculated from the angle between the incident and outgoing trajectory unit vectors:

**θ = cos⁻¹(v̂_A · v̂_B)**

The dot product is constrained to the range:

**−1 ≤ v̂_A · v̂_B ≤ 1**

to ensure numerical stability.

The resulting event-by-event scattering angles form the experimental angular distribution.

---

### 7. Angular Distribution and Gaussian Fitting

The scattering angles are histogrammed using 120 bins.

The central region of the distribution is characterized using a Gaussian function with a constant background:

**f(θ) = A exp[−½((θ − μ)/σ)²] + C**

where:

- **A** = amplitude
- **μ** = Gaussian centroid
- **σ** = Gaussian width
- **C** = constant background

The Full Width at Half Maximum (FWHM) is calculated as:

**FWHM = 2√(2 ln 2) σ ≈ 2.35482σ**

The fitted parameters are used to characterize the scattering-angle distribution at each projectile energy.

---

## Repository Structure

```text
carbon-ion-aluminum-scattering/
│
├── README.md
├── LICENSE
├── .gitignore
├── requirements.txt
│
├── data/
│   ├── raw/
│   │   ├── 21.0_MeV_per_u/
│   │   ├── 28.5_MeV_per_u/
│   │   ├── 83.0_MeV_per_u/
│   │   ├── 102.5_MeV_per_u/
│   │   ├── 107.0_MeV_per_u/
│   │   └── 110.0_MeV_per_u/
│   │
│   └── processed/
│       ├── 21.0_MeV_per_u/
│       ├── 28.5_MeV_per_u/
│       ├── 83.0_MeV_per_u/
│       ├── 102.5_MeV_per_u/
│       ├── 107.0_MeV_per_u/
│       └── 110.0_MeV_per_u/
│
├── analysis/
│   ├── 21.0_MeV_per_u.py
│   ├── 28.5_MeV_per_u.py
│   ├── 83.0_MeV_per_u.py
│   ├── 102.5_MeV_per_u.py
│   ├── 107.0_MeV_per_u.py
│   └── 110.0_MeV_per_u.py
│
├── simulations/
│   └── geant4/
│       ├── README.md
│       ├── CMakeLists.txt
│       ├── include/
│       ├── src/
│       ├── macros/
│       └── output/
│
├── figures/
│   ├── experimental/
│   │   ├── 21.0_MeV_per_u/
│   │   ├── 28.5_MeV_per_u/
│   │   ├── 83.0_MeV_per_u/
│   │   ├── 102.5_MeV_per_u/
│   │   ├── 107.0_MeV_per_u/
│   │   └── 110.0_MeV_per_u/
│   │
│   └── simulation/
│
├── results/
│   ├── scattering_angles.csv
│   ├── energy_comparison.csv
│   └── README.md
│
└── docs/
    ├── methodology.md
    └── data_description.md
```

---

## Software Requirements

The experimental data analysis was performed using Python.

Required Python packages:

```text
numpy
pandas
matplotlib
scipy
```

The Monte Carlo simulation component was developed using GEANT4.

---

## Repository Contents

- **data/raw/** — Original experimental CSV datasets.
- **data/processed/** — Processed datasets generated during the analysis.
- **analysis/** — Analysis notebooks for the six projectile energies.
- **simulations/geant4/** — GEANT4 simulation source code and related files.
- **figures/experimental/** — Figures generated from the CR-39 experimental analysis.
- **figures/simulation/** — Figures generated from the GEANT4 simulations.
- **results/** — Final scattering-angle datasets and energy-dependent comparison results.
- **docs/** — Additional documentation describing the methodology and datasets.

---

## Author

**Md. Shamsul Alam Mahfuj**

M.Sc. in Physics  
Department of Physics  
University of Chittagong, Bangladesh

---

## License

This project is distributed under the license specified in the `LICENSE` file.
