# IV_Carbon

**2D Kinetic Monte Carlo simulator for charge transport in nanostructures**  
Computes current–voltage characteristics, diffusion coefficient, mobility, and conductivity using a site‑hopping model. The electronic band gap is obtained either from a tight‑binding calculation on the supercell (automatic) or provided manually.

## Features

- **Kinetic Monte Carlo (KMC)** with Miller‑Abrahams‑like hopping rates  
  – Localization length modulated by gap, temperature, and magnetic field  
  – Source–drain geometry with user‑defined electrode thickness  
- **Tight‑binding gap** (2D)  
  – Builds Hamiltonian from structure and TB parameters  
  – Scans Brillouin zone to find minimum gap between consecutive eigenvalues  
- **Transport properties**  
  – Conductance, current, conductivity, diffusion coefficient, transit time, effective mobility  
- **Bias sweep** (automatic or single point)  
- **Flexible input** via `param.txt`  

## Requirements

- Python 3.8+
- `numpy`, `matplotlib`, `pandas`, `ase`, `pymatgen`

Install with:  
```bash
pip install numpy matplotlib pandas ase pymatgen

## Input file: `param.txt`

All simulation parameters are defined in a plain text file. Example:

# ------------------------------------------------------------
# Structure
cif_file = structure.cif
nx = 2
ny = 2
adnz = 1                 # repetitions along z (usually 1 for 2D)

# Geometry
direction = X            # transport direction (X or Y)
neighbor_cutoff = 3.5    # Å
electrode_thickness = 5.0  # Å

# Band gap
bandgap = auto           # "auto" or a number (eV)
# TB parameters (used only if bandgap = auto)
TB_T0 = -2.7             # eV
TB_D0 = 1.42             # Å
TB_BETA = 3.37
TB_TMIN = 1e-6           # eV
TB_RCUT = 4.0            # Å
TB_ONSITE = 0.0          # eV
TB_NSEG = 20             # k‑points per direction

# KMC simulation
n_balls = 1000
temperature = 300        # K
xi0 = 5.0                # Å (localization length at zero gap, reference T, B=0)
alpha_gap = 0.5          # gap dependence: xi = xi0 / (1 + alpha_gap * gap)
beta = 0.1               # temperature dependence: xi *= 1 + beta*(T/T_ref - 1)
nu0 = 1e12               # attempt frequency (Hz)
Bfield = 0.0             # Tesla
B0 = 1.0                 # reference field for magnetic suppression

# Voltage sweep
do_sweep = yes
sweep_Vmin = -1.0
sweep_Vmax = 1.0
sweep_dV = 0.05
# For single bias (if do_sweep = no):
Vbias = 0.5

## Running the code

python IV_Carbon.py

The program will:  
1. Read `param.txt`  
2. Load and replicate the CIF structure  
3. Compute the TB gap (if `bandgap = auto`)  
4. Build neighbor lists and geometry  
5. Run KMC for each bias point  
6. Save results to `sweep_vbias_results.csv`

## Output CSV columns

| Column   | Description                                      |
|----------|--------------------------------------------------|
| V        | Bias voltage (V)                                 |
| T        | Transmission coefficient (average passed balls)  |
| G        | Conductance (S)                                  |
| I        | Current (A)                                      |
| sigma    | Conductivity (S/m)                               |
| D        | Diffusion coefficient (m²/s)                     |
| tau_tr   | Mean transit time (s)                            |
| mu_eff   | Effective mobility (m²/(V·s))                    |

## How it works (brief)

- **Hopping rate** between sites `i` and `j` (separated by distance `r`):
  Γ = ν₀ * exp(-2r/ξ) / (1+2Δ) * f(ΔE)
  where ξ is the localization length, Δ the TB gap, ΔE the energy difference (from applied bias), and f the Boltzmann factor if ΔE>0.

- **KMC** injects carriers from the source electrode and tracks them until they either reach the drain or exceed a time limit.

- **TB gap** is computed on the **replicated supercell** (the actual simulation cell) to capture the electronic structure of the device region.

## Authors

J. P. Dadario Pereira, Raphael Tromer, Luiz A. Ribeiro Junior, and Douglas S. Galvao†

† Applied Physics Department, 'Gleb Wataghin' Institute of Physics, State University of Campinas, Campinas, SP, 13083-970, Brazil  
‡ University of Brasília, Institute of Physics, Brasília, Federal District, Brazil  
¶ Computational Materials Laboratory, LCCMat, Institute of Physics, University of Brasília, 70910-900, Brasília, Federal District, Brazil
