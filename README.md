# TMCM-MnCl3 Raman Mode Classifier

Classify the computed Raman-active vibrational modes of the hybrid
ferroelectric perovskite **TMCM-MnCl3**
(TMCM = trimethylchloromethylammonium, [(CH3)3NCH2Cl])
from VASP phonon data.

The TMCM-MnCl3 in Cc phase has no centrosymmetry, so all 3N−3 = 123 vibrational modes are
Raman-active. This repository assigns each mode a physical label — e.g.
_C–H stretch_, _CH3 twisting_, _CH2 wagging_, _Mn–Cl stretching_ — by analysing
the **atomic displacement pattern (eigenvector)** of every mode against the
relaxed crystal structure.

---

## Importance

A DFT phonon calculation gives you a list of modes, each with a
frequency and a Raman activity, but the mode numbers carry no physical meaning.
Before you can compare theory to an experimental Raman spectrum you need to know
**what each peak is**: i.e., is the 3 000 cm⁻¹ peak a C–H stretch? Is the 250 cm⁻¹
peak a symmetric Mn–Cl breathing mode? This tool answers that.

Thus for every mode we look at *how the atoms move* and decide, using bond/angle
metrics, which internal coordinate (stretch, bend, wag, rock, twist, umbrella, …)
dominates the motion.

---

## Repository layout

```
tmcm-mncl3-raman/
├── README.md                     <- brief intro
├── requirements.txt              <- pip dependencies
├── environment.yml               <- conda environment (alternative to pip)
├── raman_mode_classifier.py      <- the mode classifier script (main)
├── raman.ipynb                   <- same code for Jupyter notebook
├── data/                         <- example input for TMCM-MnCl3    
│   ├── POSCAR                    <- relaxed DFT structure (42 atoms, VASP5)
│   ├── vasp_raman.dat            <- mode frequencies + Raman activities
│   ├── phonon_displacements/     <- one file per mode: 123 x {mode}.txt
│   │                                  (columns: index dx dy dz)
│   └──
├── scripts/
│   ├── Extract-displacements.ipynb  <- generate phonon_displacements/ folder with the disptacement of each of the 123 modes
│   └── make_notebook.py              <- rebuild the .ipynb from the .py        
```

---


## Quick start


### 1. Run the classifier

The script reads its four inputs from the **current working directory** and
writes three output files there too.  Copy `raman_mode_classifier.py` (or the
notebook) into the folder that already contains your VASP files, then run:

```bash
cd /path/to/your/TMCM-MnCl3/DFT        # must contain POSCAR, vasp_raman.dat, etc.
cp /path/to/repo/raman_mode_classifier.py .
python raman_mode_classifier.py
```

If you want symmetry labels (`A'` / `A''`), place `OUTCAR.phon` in the same
directory (or symlink it):

```bash
ln -s /path/to/your/OUTCAR.phon .
```

> **Quick test inside the repo.** The `data/` directory holds sample copies of
> the inputs.  To run without a 2.4 GB OUTCAR symlinked the A'/A'' column
> simply shows `?`; with it:
>
> ```bash
> cd data
> python ../raman_mode_classifier.py
> ```

### 3. What comes out

| File | Content |
|------|---------|
| `phonon_mode_report.txt` | human-readable, every mode + its top assignments |
| `phonon_mode_summary.csv` | one row per mode: number, cm⁻¹, activity, symmetry, 4 assignments |
| `raman_table.tex` | LaTeX table for your paper |

Sample of `phonon_mode_summary.csv`:

```
Mode,Frequency_cm-1,Raman_Activity,Symmetry,Assign_1,Assign_2,...
15,2999.91,6455.91,$A'$,"CH3 symmetric stretching (C5)",...
68,788.84,392.60,$A''$,"C-Cl asymmetric stretching (C1-Cl7)",...
```

### 4. Inspect the results

```bash
head -20 outputs/phonon_mode_report.txt
```

or open the notebook:

```bash
jupyter notebook raman.ipynb
```

---

## Input data formats

### `data/vasp_raman.dat`

Produced by the phonopy [`vasp_raman.py`](https://phonopy.github.io/phonopy/) accessory (Raman-tensor post-processing).
Five columns (`#` lines are comments):

```
# mode    freq(cm-1)    alpha    beta2    activity
124     0.22823   0.0000000   0.0000016   0.0000114
...
```

* `alpha`  — trace of the Raman tensor (isotropic polarizability derivative), Å²/amu-style units
* `beta2 = beta²` — anisotropy squared of the Raman tensor
* `activity = 45·alpha² + 7·beta²` — the powder (orientationally averaged) Raman scattering factor

### `data/phonon_displacements/{mode}.txt`

Numeric columns per line:

```
# mode 90
# index dx dy dz
1 -0.001090 0.000436 -0.000493
...
```

Each line is one atom (`index` = 1..42, matching the POSCAR atom order),
and `dx dy dz` is its **displacement vector** (Å) in that mode. These come
straight from the *"Eigenvectors and eigenvalues of the dynamical matrix"*
section of `OUTCAR.phon`. Regenerate them with:

```bash
python scripts/extract_modes_from_outcar.py /path/to/OUTCAR.phon data/phonon_displacements
```

> The extracted files are already committed so you can run the pipeline
> without owning a 2.4 GB OUTCAR. Only the symmetry column needs the real
> OUTCAR.

### `data/POSCAR`

VASP5 structure file (42 atoms: 2 Mn + 8 Cl + 2 N + 8 C + 22 H). The species
must be on line 6 (`Mn Cl N C H`) exactly as in this repo.

---

## How the classification works (short version)

1. **Read the structure.** `load_poscar()` parses the lattice + coordinates,
   wraps all positions into the unit cell (periodic boundary conditions).

2. **Build the bond graph.** Two atoms are neighbours if their
   minimum-image distance is below a per-element-pair cutoff
   (`BOND_CUTOFFS`). This gives Mn a 6-fold Cl octahedron, the N four C
   neighbours (TMCM = N(CH3)3·CH2Cl), the two CH2Cl carbons an
   N+2H+Cl neighbourhood, and the six CH3 carbons an N+3H neighbourhood.

3. **For each mode**:
   - Load its displacement vector (`dx dy dz` per atom).
   - Measure **bond-length changes** (`Δr`) between every bonded pair.
   - Measure **angle changes** (`Δθ`) at every relevant centre
     (H-C-H, Cl-C-N, Cl-Mn-Cl, …).
   - Project H-pair motion onto a local frame (`n`, `t`) at the carbon
     to separate **wagging** (out-of-plane, in phase), **twisting**
     (out-of-plane, antiphase), **rocking** (in-plane, in phase) from
     **scissoring** (in-plane, antiphase, caught via the H-C-H angle).

4. **Thresholds decide what counts** (see table below).

5. **Rank and report.** Assignments are sorted by an internal *magnitude*
   and the top few written to the CSV/LaTeX outputs.

### The thresholds

| Constant | Value | Meaning | What it implies |
|----------|-------|---------|-----------------|
| `FREQ_MIN` / `FREQ_MAX` | 5 / 3200 cm⁻¹ | only analyse this frequency window | drops the 3 acoustic modes (0.2–0.5 cm⁻¹) and nothing else |
| `CH_STRETCH_MIN` | 2800 cm⁻¹ | C–H stretch is only labelled above this | prevents a bend/skeleton mode being called a stretch |
| `DISP_ADAPTIVE_FRACTION` | 0.02 | adaptive floor = 2 % of the *largest atomic displacement* in the mode | secondary participants need to move ≥ 2 % of the strongest mover |
| `BOND_ABS_THRESHOLD` | 0.01 Å | absolute bond-change floor for a "stretch" label | a bond that stretches less than 10⁻² Å is ignored |
| `ANGLE_CHANGE_THRESHOLD` | 10° | angle-change floor for bend/scissor/wag labels | angles that move < 10° are ignored |
| `RELATIVE_ASSIGNMENT_THRESHOLD` | 0.30 | a LaTeX assignment needs magnitude ≥ 30 % of the mode's top magnitude | keeps the table meaningful, discards weak caveats |
| `MAX_LATEX_ASSIGNMENTS` | 4 | max labels per LaTeX row | keeps the table compact |
| `BOND_CUTOFFS` | C–H 1.1, N–C 1.52, Mn–Cl 2.6, C–Cl 1.78, Cl–Cl 3.13 Å | neighbour definition | tuned to covalent radii; defines the bond graph |

> **Units note (fixed).** Length-based motions (bond stretches, wag, rock, twist,
> umbrella) are ranked in Angstrom; angle changes (scissoring, bend) are in
> degrees. The ranking keeps these two sets separate so that a 20 degscissor no
> longer crowds out a 0.5 A wag.  In the LaTeX table the 0.30 relative threshold
> is applied only within the length set; angle-based assignments then fill the
> remaining slots up to `MAX_LATEX_ASSIGNMENTS`.

---

## What the assignment labels mean
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

