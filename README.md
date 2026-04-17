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
