# TMCM-MnCl3 Raman Mode Classifier

Classify the computed Raman-active vibrational modes of the hybrid
ferroelectric perovskite **TMCM-MnCl3**
(TMCM = trimethylchloromethylammonium, [(CH3)3NCH2Cl])
from VASP phonon data.

The TMCM-MnCl3 in Cc phase has no centrosymmetry, so all 3N−3 = 123 vibrational modes are
Raman-active. This repository assigns each mode a physical label — e.g.
_C–H stretch_, _CH3 twisting_, _CH2 wagging_, _Mn–Cl stretching_ — by analysing
the eigenvector of every mode against the relaxed POSCAR.

---

## Importance

A DFT phonon calculation gives you a list of modes, each with a
frequency and a Raman activity, but the mode numbers carry no physical meaning.
Before you can compare theory to an experimental Raman spectrum you need to know
**what each peak is**: i.e., is the 3 000 cm⁻¹ peak a C–H stretch? Is the 250 cm⁻¹
peak a symmetric Mn–Cl breathing mode? This tool answers that.

Thus for every mode we look at *how the atoms move* and decide, using bond/angle
metrics, what (stretch, bend, wag, rock, twist...) dominates the motion.

---

## Repository layout

```
tmcm-mncl3-raman/
├── README.md                     <- brief intro
├── requirements.txt             
├── data/                         <- example input for TMCM-MnCl3    
│   ├── POSCAR                    <- relaxed DFT structure (42 atoms, VASP5)
│   ├── vasp_raman.dat            <- mode frequencies + Raman activities
│   ├── phonon_displacements/     <- generate from Extract-displacements.ipynb                         
│   └──
├── scripts/
│   ├── Extract-displacements.ipynb  <- generate phonon_displacements/ folder with the disptacement of each of the 123 modes
│   └── raman.ipynb      <- the mode classifier Jupyter notebook script (main)


```

---

## Input data formats

### `data/vasp_raman.dat`

Produced using from the DFT OUTCAR using https://github.com/raman-sc/VASP.
Five columns :

```
# mode    freq(cm-1)    alpha    beta2    activity
...
```

* `alpha`  — trace of the Raman tensor (isotropic polarizability derivative), Å²/amu-style units
* `beta2 = beta²` — anisotropy squared of the Raman tensor
* `activity = 45·alpha² + 7·beta²` — the powder (orientationally averaged) Raman scattering factor

### `data/phonon_displacements/{mode}.txt`

Numeric columns per line:

```

# atom_index dx dy dz
...
```

Each line is the atom index matching the POSCAR atom order and
`dx dy dz` is its **displacement vector** (Å) in that mode. These come
straight from the *"Eigenvectors and eigenvalues of the dynamical matrix"*
section of `OUTCAR`. Regenerate them with scripts/Extract-displacements.ipynb



## How the classification works

1. **Read the structure.** `load_poscar()` parses the lattice + coordinates,
   wraps all positions into the unit cell (periodic boundary conditions).

2. **Build the bond graph.** Two atoms are neighbours if their
   minimum-image distance is below a cutoff
   (`BOND_CUTOFFS`). This gives Mn a 6-Cl octahedra, the N four C
   neighbours, the two CH2Cl carbons an
   N+2H+Cl neighbourhood, and the six CH3 carbons an N+3H neighbourhood.

3. **For each mode**:
   - Load its displacement vector (`dx dy dz` per atom).
   - Measure **bond-length changes** (`Δr`) between every bonded pair.
   - Measure **angle changes** (`Δθ`) at every relevant centre
     (H-C-H, Cl-C-N, Cl-Mn-Cl, …).

4. **Thresholds**.
   BOND_CUTOFFS - tuned to covalent radii; defines the bond graph
   ANGLE_CHANGE_THRESHOLD - minimum angle-change for bend/scissor/wag labels. Angles that move < 10° are ignored
   BOND_ABS_THRESHOLD -minimum absolute bond-change for a "stretch" label. A bond that stretches less than 10⁻² Å is ignored |
   
5. **Rank modes.** Assignments are sorted by an internal *magnitude*
   and the top few written to the CSV/LaTeX outputs.


---

## Assignment labels Legend
ν - stretching
ω - wagging
τ - twisting
ρ - rocking
β - Scissoring
as - asymmetric
s - symmetric

## Reproducing everything from scratch

```bash
# 1. extract per-mode displacement files from your OUTCAR from DFPT calculations
- Extract-displacements.ipynb 

# 2. classify the modes
python raman_mode_classifier.py

```

---

## Requirements / version notes

- Python ≥ 3.8
- `numpy` (required)
- `scipy`, `matplotlib` (required only for plotting in the notebook / docs)
- The `data/` files in this repo were produced with VASP 6.3.1 + phonopy.

---

## License

