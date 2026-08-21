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

The datasets contain the centroid coordinates of etched ion tracks extracted from PitFit. The primary coordinates used in the analysis are:

- `GravX`
- `GravY`

The coordinates are expressed in micrometers (µm).

---

# Analysis Workflow

The same analysis procedure was applied independently to all six beam energies.

## 1. Data Loading

For each energy, four CSV files corresponding to the four detector surfaces are loaded:

```text
A-front
A-back
B-front
B-back
```

Entries with missing values in `GravX` or `GravY` are removed before further analysis.

---

## 2. Front-Back Alignment

The front and back surfaces of a CR-39 detector may have different coordinate origins due to scanning and measurement offsets.

For each detector sheet, coordinate differences are calculated as:

\[
\Delta x = x_{\mathrm{back}} - x_{\mathrm{front}}
\]

\[
\Delta y = y_{\mathrm{back}} - y_{\mathrm{front}}
\]

A two-dimensional Kernel Density Estimation (KDE) is applied to the \((\Delta x,\Delta y)\) distribution.

The peak of the KDE distribution gives the relative translational offset between the front and back coordinate systems.

The back-surface coordinates are then corrected using:

\[
x_{\mathrm{back}}' =
x_{\mathrm{back}} -
\Delta x_{\mathrm{peak}}
\]

\[
y_{\mathrm{back}}' =
y_{\mathrm{back}} -
\Delta y_{\mathrm{peak}}
\]

This alignment procedure is performed independently for both Sheet A and Sheet B.

---

## 3. Front-Back Track Matching

After alignment, tracks on the front and back surfaces of each CR-39 sheet are matched.

A KD-tree nearest-neighbour search is used to identify candidate tracks.

The matching radius is determined from the spatial spread of the alignment distribution:

\[
r = 4\sigma
\]

Only unique one-to-one matches are retained. Ambiguous matches and tracks without a valid counterpart are excluded.

This produces reconstructed track pairs for:

- Sheet A: incident trajectories before the aluminum target.
- Sheet B: outgoing trajectories after the aluminum target.

---

## 4. Three-Dimensional Trajectory Reconstruction

For every matched front-back track pair, a three-dimensional trajectory vector is reconstructed.

For Sheet A:

\[
\vec{v}_A =
\left(
x_{A,b}-x_{A,f},
y_{A,b}-y_{A,f},
t_A
\right)
\]

For Sheet B:

\[
\vec{v}_B =
\left(
x_{B,b}-x_{B,f},
y_{B,b}-y_{B,f},
t_B
\right)
\]

where \(t_A\) and \(t_B\) are the post-etch thicknesses of the corresponding CR-39 sheets.

The unit vectors are calculated as:

\[
\hat{v} =
\frac{\vec{v}}
{|\vec{v}|}
\]

The trajectory reconstructed in Sheet A represents the direction of the carbon ion before passing through the aluminum target, while the trajectory reconstructed in Sheet B represents the direction after transmission through the target.

---

## 5. Cross-Matching of Incident and Outgoing Tracks

The trajectory reconstructed in Sheet A is projected to the plane of the front surface of Sheet B.

The projected position is calculated using the reconstructed direction and the known separation between the detector sheets.

A KD-tree nearest-neighbour search is then used to identify the corresponding trajectory in Sheet B.

Only successfully cross-matched trajectories are retained for the final scattering-angle calculation.

---

## 6. Scattering-Angle Calculation

For each cross-matched event, the scattering angle is calculated from the angle between the incident and outgoing trajectory vectors:

\[
\theta =
\cos^{-1}
\left(
\hat{v}_A
\cdot
\hat{v}_B
\right)
\]

The dot product is constrained to the range:

\[
[-1,1]
\]

to ensure numerical stability.

The resulting event-by-event scattering angles form the experimental angular distribution.

---

## 7. Angular Distribution and Gaussian Fitting

The scattering angles are histogrammed using 120 bins.

The central region of the distribution is characterized using a Gaussian function with a constant background:

\[
f(\theta) =
A
\exp
\left[
-\frac{1}{2}
\left(
\frac{\theta-\mu}{\sigma}
\right)^2
\right]
+
C
\]

where:

- \(A\) is the amplitude,
- \(\mu\) is the Gaussian centroid,
- \(\sigma\) is the Gaussian width,
- \(C\) is a constant background.

The full width at half maximum is calculated as:

\[
\mathrm{FWHM}
=
2\sqrt{2\ln2}\sigma
\approx
2.35482\sigma
\]

The fitted parameters are used to characterize the scattering-angle distribution at each projectile energy.

---

# Repository Structure

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
│   ├── 21.0_MeV_per_u.ipynb
│   ├── 28.5_MeV_per_u.ipynb
│   ├── 83.0_MeV_per_u.ipynb
│   ├── 102.5_MeV_per_u.ipynb
│   ├── 107.0_MeV_per_u.ipynb
│   └── 110.0_MeV_per_u.ipynb
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

# Software Requirements

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

# Contents

- **data/raw/** contains the original experimental CSV datasets.
- **data/processed/** contains processed datasets generated during the analysis.
- **analysis/** contains the analysis notebooks for each projectile energy.
- **simulations/geant4/** contains the GEANT4 simulation source code and related files.
- **figures/experimental/** contains figures generated from the CR-39 experimental analysis.
- **figures/simulation/** contains figures generated from GEANT4 simulations.
- **results/** contains the final scattering-angle datasets and energy-dependent comparison results.
- **docs/** contains additional documentation describing the methodology and datasets.

---

# Author

**Md. Shamsul Alam Mahfuj**

M.Sc. in Physics  
Department of Physics  
University of Chittagong, Bangladesh

---

# License

This project is distributed under the license specified in the `LICENSE` file.
