IV_Carbon
<p align="center"> <b>2D Kinetic Monte Carlo simulator for charge transport in nanostructures</b><br> Site-hopping transport • Tight-binding gap • Full I–V characterization </p> <p align="center"> <img src="https://img.shields.io/badge/python-3.8%2B-blue.svg"> <img src="https://img.shields.io/badge/status-active-success.svg"> <img src="https://img.shields.io/badge/license-MIT-green.svg"> <img src="https://img.shields.io/badge/domain-nanotransport-orange.svg"> </p>
Overview

IV_Carbon computes charge transport properties in 2D nanostructures using a Kinetic Monte Carlo (KMC) framework with Miller–Abrahams-type hopping.

It supports:

Automatic tight-binding band gap extraction
Full I–V curve simulation
Transport observables such as mobility, diffusion, and conductivity
Features
Transport Model
Kinetic Monte Carlo (KMC)
Miller–Abrahams-like hopping rates
Localization length dependent on:
band gap
temperature
magnetic field
Electronic Structure
Tight-binding Hamiltonian (2D)
Brillouin zone sampling
Automatic gap detection
Device Physics
Source–drain geometry
User-defined electrode thickness
Voltage sweep or single bias
Outputs
Current–voltage (I–V)
Conductance and conductivity
Diffusion coefficient
Transit time
Effective mobility
Installation
Requirements
Python 3.8+
numpy
matplotlib
pandas
ase
pymatgen
Install dependencies
pip install numpy matplotlib pandas ase pymatgen
Input File: param.txt

All simulation parameters are defined here.

<details> <summary><b>Example configuration</b></summary>
# ------------------------------------------------------------
# Structure
cif_file = structure.cif
nx = 2
ny = 2
adnz = 1

# Geometry
direction = X
neighbor_cutoff = 3.5
electrode_thickness = 5.0

# Band gap
bandgap = auto

# TB parameters
TB_T0 = -2.7
TB_D0 = 1.42
TB_BETA = 3.37
TB_TMIN = 1e-6
TB_RCUT = 4.0
TB_ONSITE = 0.0
TB_NSEG = 20

# KMC
n_balls = 1000
temperature = 300
xi0 = 5.0
alpha_gap = 0.5
beta = 0.1
nu0 = 1e12
Bfield = 0.0
B0 = 1.0

# Voltage sweep
do_sweep = yes
sweep_Vmin = -1.0
sweep_Vmax = 1.0
sweep_dV = 0.05

# Single bias (if do_sweep = no)
Vbias = 0.5
</details>
Running the Simulation
python IV_Carbon.py
Workflow
param.txt → CIF structure → Supercell → TB gap → Neighbor graph → KMC → CSV output

Steps:

Read input parameters
Load and replicate structure
Compute TB band gap (if needed)
Build neighbor network
Run KMC transport
Save results
Output
File
sweep_vbias_results.csv
Columns
Variable	Meaning
V	Bias voltage (V)
T	Transmission coefficient
G	Conductance (S)
I	Current (A)
sigma	Conductivity (S/m)
D	Diffusion coefficient (m²/s)
tau_tr	Transit time (s)
mu_eff	Mobility (m²/V·s)
Model Details
Hopping Rate
Γ=ν
0
	​

⋅e
−2r/ξ
⋅
1+2Δ
1
	​

⋅f(ΔE)
r: distance between sites
ξ: localization length
Δ: band gap
ν
0
	​

: attempt frequency
f(ΔE): Boltzmann factor
Kinetic Monte Carlo
Carriers injected at source electrode
Perform stochastic hopping
Simulation ends when:
carrier reaches drain, or
time limit is exceeded
Tight-Binding Gap
Computed on the replicated supercell
Hamiltonian built from atomic positions + TB parameters
Gap = minimum difference between consecutive eigenvalues across k-space
Project Structure
IV_Carbon/
├── IV_Carbon.py
├── param.txt
├── structure.cif
└── results/
Roadmap
 Parallel KMC (MPI / multiprocessing)
 Disorder models (energetic / positional)
 GUI for param setup
 3D extension
Authors
J. P. Dadario Pereira
Raphael Tromer
Luiz A. Ribeiro Junior
Douglas S. Galvão

Affiliations:

UNICAMP – Gleb Wataghin Institute of Physics
University of Brasília (UnB)
LCCMat – Computational Materials Laboratory
Citation

If you use this code in research, consider citing:

@software{iv_carbon,
  title = {IV\_Carbon: KMC Transport Simulator},
  author = {Dadario Pereira, J. P. and collaborators},
  year = {2026}
}
