# 2D Galton atomic simulation: KMC + IV + diffusion + mobility + TB gap
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from pymatgen.core import Structure
from dataclasses import dataclass
import pandas as pd
import warnings

warnings.filterwarnings("ignore")

# Physical constants
ELECTRON_CHARGE = 1.602176634e-19   # C
KB_EV = 8.617333262e-5             # eV/K
KB_J = 1.380649e-23                # J/K
PLANCK_J = 6.62607015e-34          # J·s

# Numerical limits
MIN_GAP = 1e-6      # eV
SWEEP_EPS = 1e-12
TEMP_REF = 300.0    # K

# ----------------------------------------------------------------------
# Parameter loading
# ----------------------------------------------------------------------
def load_params(fname="param.txt"):
    params = {}
    with open(fname, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, val = line.split("=", 1)
            key = key.strip()
            val = val.split("#", 1)[0].strip().replace('"', "").replace("'", "").replace(",", ".")
            if val.lower() in ("yes", "no"):
                params[key] = val.lower()
            elif val.upper() in ("X", "Y"):
                params[key] = val.upper()
            else:
                try:
                    params[key] = float(val) if ("." in val or "e" in val.lower()) else int(val)
                except ValueError:
                    params[key] = val
    return params

# ----------------------------------------------------------------------
# Geometry utilities
# ----------------------------------------------------------------------
def read_and_repeat(cif_file, nx, ny, nz):
    return read(cif_file).repeat((nx, ny, nz))

def get_positions(atoms):
    pos = atoms.get_positions()
    return pos[:, 0], pos[:, 1]

def align_to_x(X, Y, cell):
    a1 = cell[0]
    angle = np.arctan2(a1[1], a1[0])
    rot = np.array([[np.cos(-angle), -np.sin(-angle)],
                    [np.sin(-angle),  np.cos(-angle)]])
    xy = np.vstack([X, Y]).T @ rot.T
    return xy[:, 0], xy[:, 1]

def get_dimensions(X, Y, direction):
    if direction == "X":
        return X.max() - X.min(), Y.max() - Y.min()
    else:
        return Y.max() - Y.min(), X.max() - X.min()

# ----------------------------------------------------------------------
# Neighbors
# ----------------------------------------------------------------------
def find_neighbors(X, Y, cutoff):
    xy = np.vstack([X, Y]).T
    nbrs = [[] for _ in range(len(xy))]
    for i in range(len(xy)):
        for j in range(i+1, len(xy)):
            d = np.hypot(*(xy[j] - xy[i]))
            if d <= cutoff:
                nbrs[i].append((j, d))
                nbrs[j].append((i, d))
    return nbrs

# ----------------------------------------------------------------------
# Site energies under bias
# ----------------------------------------------------------------------
def get_site_energies(X, Y, voltage, direction):
    coord = X if direction == "X" else Y
    span = coord.max() - coord.min()
    if span == 0:
        return np.zeros_like(coord)
    return -voltage * (coord - coord.min()) / span

# ----------------------------------------------------------------------
# Hopping rates (KMC)
# ----------------------------------------------------------------------
def calc_hop_rates(src_idx, nbr_list, E_src, E_dst, gap, T, xi0, alpha, beta, nu0, B, B0):
    xi = xi0 / (1 + alpha * gap)
    xi *= 1 + beta * (T / TEMP_REF - 1)
    xi /= np.sqrt(1 + (B / B0)**2)

    targets = []
    rates = []
    for j, dist in nbr_list:
        dE = E_dst[j] - E_src[src_idx]
        base = nu0 * np.exp(-2 * dist / xi) / (1 + 2 * gap)
        f = np.exp(-dE / (KB_EV * T)) if (dE > 0 and T > 0) else 1.0
        rate = base * f
        if rate > 0:
            targets.append(j)
            rates.append(rate)
    return np.array(targets, dtype=int), np.array(rates, dtype=float)

# ----------------------------------------------------------------------
# KMC simulation
# ----------------------------------------------------------------------
def run_kmc(X, Y, gap, voltage, direction, nballs, nbrs, T, p, elec_thick, tmax=1e3):
    E = get_site_energies(X, Y, voltage, direction)
    # local copies to avoid aliasing
    E_src = E
    E_dst = E.copy()

    if direction == "X":
        source_sites = np.where(X < X.min() + elec_thick)[0]
        coord = X
    else:
        source_sites = np.where(Y < Y.min() + elec_thick)[0]
        coord = Y

    drain_thresh = coord.max() - elec_thick

    results = []
    for _ in range(nballs):
        site = np.random.choice(source_sites)
        x0 = coord[site]
        t = 0.0
        passed = False

        while t < tmax:
            targets, rates = calc_hop_rates(
                site, nbrs[site], E_src, E_dst, gap, T,
                p["xi0"], p["alpha_gap"], p["beta"], p["nu0"],
                p["Bfield"], p["B0"]
            )
            if len(rates) == 0:
                break

            W = rates.sum()
            dt = -np.log(np.random.rand()) / W
            t += dt

            rnd = np.random.rand() * W
            chosen = np.searchsorted(np.cumsum(rates), rnd)
            chosen = min(chosen, len(targets)-1)
            site = targets[chosen]

            if coord[site] > drain_thresh:
                passed = True
                break

        dx = coord[site] - x0
        results.append((passed, t, dx))

    return results

# ----------------------------------------------------------------------
# Transport properties
# ----------------------------------------------------------------------
def compute_transport(passed, voltage, length, area):
    G0 = 2 * ELECTRON_CHARGE**2 / PLANCK_J
    trans = float(np.mean(passed)) if len(passed) else 0.0
    G = G0 * trans
    I = G * voltage
    sigma = G * area / length if length > 0 else 0.0
    return trans, G, I, sigma

# ----------------------------------------------------------------------
# Tight‑binding gap (2D supercell, consecutive eigenvalues)
# ----------------------------------------------------------------------
class TBModel:
    def __init__(self, t0, d0, beta, tmin):
        self.t0 = t0
        self.d0 = d0
        self.beta = beta
        self.tmin = tmin
    def hopping(self, dr):
        val = self.t0 * np.exp(-self.beta * (np.linalg.norm(dr) / self.d0 - 1))
        return 0.0 if abs(val) < self.tmin else float(val)

def get_tb_pairs(structure, rcut):
    pairs = []
    for i, site in enumerate(structure):
        for nbr in structure.get_neighbors(site, rcut, include_index=True):
            if nbr.index > i:
                dr = np.array(nbr.coords) - np.array(site.coords)
                pairs.append((i, nbr.index, dr))
    return pairs

def build_tb_h(structure, pairs, model, onsite, kpt):
    n = len(structure)
    H = np.eye(n, dtype=complex) * onsite
    kcart = structure.lattice.reciprocal_lattice.get_cartesian_coords(kpt)
    for i, j, dr in pairs:
        t = model.hopping(dr)
        if t == 0.0:
            continue
        phase = np.exp(1j * np.dot(kcart, dr))
        H[i, j] += t * phase
        H[j, i] += np.conj(t * phase)
    return H

def compute_tb_gap(structure, p):
    model = TBModel(p["TB_T0"], p["TB_D0"], p["TB_BETA"], p["TB_TMIN"])
    pairs = get_tb_pairs(structure, p["TB_RCUT"])
    ks = np.linspace(0, 1, int(p["TB_NSEG"]), endpoint=False)
    min_gap = np.inf
    for kx in ks:
        for ky in ks:
            ev = np.sort(np.real(np.linalg.eigvalsh(
                build_tb_h(structure, pairs, model, p["TB_ONSITE"], [kx, ky, 0])
            )))
            gaps = np.diff(ev)
            pos_gaps = gaps[gaps > MIN_GAP]
            if len(pos_gaps):
                min_gap = min(min_gap, pos_gaps.min())
    return 0.0 if not np.isfinite(min_gap) else float(min_gap)

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    p = load_params("param.txt")

    # Inputs
    cif_file = p["cif_file"]
    nx, ny, nz = int(p["nx"]), int(p["ny"]), int(p["adnz"])
    direction = p["direction"]
    nballs = int(p["n_balls"])
    T = float(p["temperature"])
    elec_thick = float(p["electrode_thickness"])

    # Read structure
    atoms = read_and_repeat(cif_file, nx, ny, nz)

    # TB gap (auto or fixed)
    if isinstance(p["bandgap"], str) and p["bandgap"].strip().lower() == "auto":
        print("Calculating TB gap on supercell (consecutive ΔE)...")
        bandgap = compute_tb_gap(Structure.from_ase_atoms(atoms), p)
        print(f"TB gap = {bandgap:.6e} eV")
    else:
        bandgap = float(p["bandgap"])

    # Prepare coordinates
    X, Y = get_positions(atoms)
    X, Y = align_to_x(X, Y, atoms.get_cell())
    nbrs = find_neighbors(X, Y, float(p["neighbor_cutoff"]))
    length, area = get_dimensions(X, Y, direction)

    # Voltage sweep
    do_sweep = str(p.get("do_sweep", "yes")).lower()
    if do_sweep == "yes":
        voltages = np.arange(float(p["sweep_Vmin"]), float(p["sweep_Vmax"]) + SWEEP_EPS, float(p["sweep_dV"]))
    else:
        voltages = np.array([float(p["Vbias"])])

    rows = []
    print("Running IV sweep...")
    for V in voltages:
        sim = run_kmc(X, Y, bandgap, V, direction, nballs, nbrs, T, p, elec_thick)
        passed = np.array([r[0] for r in sim], dtype=bool)
        times = np.array([r[1] for r in sim], dtype=float)
        dx = np.array([r[2] for r in sim], dtype=float)

        trans, G, I, sigma = compute_transport(passed, V, length, area)
        t_mean = float(np.mean(times)) if np.mean(times) > 0 else np.nan
        msd = float(np.mean(dx * dx))
        D = msd / (2 * t_mean) if (t_mean > 0 and np.isfinite(t_mean)) else 0.0
        tau = t_mean if np.isfinite(t_mean) else 0.0
        mu = (ELECTRON_CHARGE * D) / (KB_J * T) if T > 0 else 0.0

        rows.append([V, trans, G, I, sigma, D, tau, mu])

    df = pd.DataFrame(rows, columns=["V", "T", "G", "I", "sigma", "D", "tau_tr", "mu_eff"])
    df.to_csv("sweep_vbias_results.csv", index=False)
    print("IV sweep completed. Results saved to sweep_vbias_results.csv")