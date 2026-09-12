#!/usr/bin/env python3
"""
TMCM-MnCl3 Phonon Mode Classification & Raman Table Generator
==============================================================

Classifies each real-frequency phonon mode of TMCM-MnCl3 into local
motions (C-H/N-H stretches, CH2/CH3 bends, C-Cl / Mn-Cl / C-N stretches,
MnCl6 deformations, ...) and writes a text report, a CSV summary, and a
LaTeX Raman table.

Input files (read from the current working directory):
  - POSCAR                  atomic structure (species line must be present)
  - vasp_raman.dat          mode index, frequency (cm^-1), Raman activity
  - OUTCAR.phon             (optional) phonon eigenvectors, used only for the
                            A'/A'' symmetry label
  - phonon_displacements/   per-mode displacement files "N.txt", one line per
                            atom in the format "index dx dy dz" (4 columns) or
                            "X Y Z dx dy dz" (6 columns)

Corrections relative to the original notebook:
  * Mixing units in the magnitude ranking has been fixed. Angle changes are
    stored in degrees and bond displacements in Angstrom; numerically these
    are not comparable (a 20 degate scissor outranks a 0.5 A wag). Results
    are now tagged by unit, length-based motions rank first, and angle-based
    motions fill the remaining top-N / LaTeX slots in order of
    |angle change|.  Previously some modes (e.g. 24-26, 33-36, 43-44) were
    rendered with only "beta(CH2)" because the scissoring angle numerically
    dominated all length terms.
  * OUTCAR.phon is parsed streaming (line by line) instead of readlines() so
    the 2.4 GB file is not loaded wholesale into memory.  Outputs of the
    parser are unchanged.
  * No mass weighting is applied: VASP eigenvectors already diagonalise the
    (mass-weighted) dynamical matrix, so raw components are used for
    ranking, as in the original script.

Known limitation (kept by design, matching the original):
  * The A'/A'' symmetry column is a heuristic (projection of the mode onto
    the lattice b vector, threshold 0.5), not a group-theory assignment.
"""

import os
import re
import sys
import csv
import numpy as np

# ============================================================
# CONFIG
# ============================================================

POSCAR_FILE = "POSCAR"
RAMAN_FILE = "vasp_raman.dat"
OUTCAR_PHON_FILE = "OUTCAR.phon"     # optional, for symmetry labels only
DISPLACEMENTS_DIR = "phonon_displacements"

OUTPUT_REPORT_FILE = "phonon_mode_report.txt"
OUTPUT_SUMMARY_CSV = "phonon_mode_summary.csv"
OUTPUT_LATEX_FILE = "raman_table.tex"

# Fallback atom species/counts, used ONLY if POSCAR has no species line (VASP4)
ATOM_TYPES_FALLBACK = ['Mn', 'Cl', 'N', 'C', 'H']
ATOM_COUNTS_FALLBACK = [2, 8, 2, 8, 22]

FREQ_MIN = 5.0
FREQ_MAX = 3200.0
CH_STRETCH_MIN = 2800.0

BOND_CUTOFFS = {
    ('C', 'H'): 1.1,
    ('N', 'C'): 1.52,
    ('Mn', 'Cl'): 2.6,
    ('C', 'Cl'): 1.78,
    ('Cl', 'Cl'): 3.13,
}

ANGLE_CHANGE_THRESHOLD = 10.0    # degrees
BOND_ABS_THRESHOLD = 1e-2        # Angstrom
DISP_ADAPTIVE_FRACTION = 0.02
RELATIVE_ASSIGNMENT_THRESHOLD = 0.30
MAX_LATEX_ASSIGNMENTS = 4
CONSOLE_TOP_N = 4
CSV_TOP_N = 4
REPORT_TOP_N = 15


def bond_cutoff(label_i, label_j):
    """Single-direction lookup for pair cutoffs."""
    return BOND_CUTOFFS.get((label_i, label_j), BOND_CUTOFFS.get((label_j, label_i)))


# ============================================================
# POSCAR LOADING (with PBC wrap + species-line parsing)
# ============================================================

def load_poscar(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"POSCAR file not found: {path}")

    with open(path) as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    if len(lines) < 7:
        raise ValueError(f"POSCAR file appears incomplete: {path}")

    scale = float(lines[1])
    if scale < 0:
        raise ValueError("Negative POSCAR scale factor is not supported.")

    lattice_raw = np.array(
        [[float(x) for x in lines[i].split()] for i in range(2, 5)]
    )
    lattice = lattice_raw * scale

    species_line = lines[5].split()
    counts_line = lines[6].split()

    try:
        counts = [int(x) for x in counts_line]
        natoms = sum(counts)
    except ValueError:
        counts = [int(x) for x in species_line]
        natoms = sum(counts)
        species_line = None

    element_labels = None
    if species_line is not None and len(species_line) == len(counts):
        element_labels = []
        for sym, cnt in zip(species_line, counts):
            element_labels += [sym] * cnt

    coord_type_line = lines[7] if species_line is not None else lines[6]
    coord_start = 8 if species_line is not None else 7
    is_direct = coord_type_line.lower().startswith('d')

    coords = []
    for i in range(coord_start, coord_start + natoms):
        parts = lines[i].split()
        coords.append([float(x) for x in parts[:3]])
    coords = np.array(coords)

    if is_direct:
        coords = coords @ lattice
    else:
        coords = coords * scale

    inv_lat = np.linalg.inv(lattice)
    frac = coords @ inv_lat
    frac = frac % 1.0
    coords = frac @ lattice

    return coords, lattice, element_labels, natoms


try:
    coords, lattice, element_labels_from_poscar, num_atoms = load_poscar(POSCAR_FILE)
except Exception as e:
    print(f"ERROR loading POSCAR: {e}", file=sys.stderr)
    sys.exit(1)

inv_lattice = np.linalg.inv(lattice)

# ============================================================
# ATOM LABELS
# ============================================================

if element_labels_from_poscar is not None:
    atom_labels = np.array(element_labels_from_poscar)
    source = "POSCAR species line"
else:
    fallback_labels = []
    for atype, count in zip(ATOM_TYPES_FALLBACK, ATOM_COUNTS_FALLBACK):
        fallback_labels += [atype] * count
    if len(fallback_labels) != num_atoms:
        print(
            f"ERROR: Fallback config sum ({len(fallback_labels)}) "
            f"does not match POSCAR atom count ({num_atoms}).",
            file=sys.stderr,
        )
        sys.exit(1)
    atom_labels = np.array(fallback_labels)
    source = "ATOM_TYPES_FALLBACK config"

print(f"Loaded {num_atoms} atoms from '{POSCAR_FILE}' (species from {source}).")

group_indices = {
    'TMCM': np.where(np.isin(atom_labels, ['C', 'N', 'H']))[0],
    'MnCl6': np.where(np.isin(atom_labels, ['Mn', 'Cl']))[0],
}
for atype in set(atom_labels):
    group_indices[atype] = np.where(atom_labels == atype)[0]

relative_indices = {}
for atype in sorted(set(atom_labels)):
    idxs = np.where(atom_labels == atype)[0]
    for ri, ai in enumerate(idxs, 1):
        relative_indices[ai] = (atype, ri)


def fmt(idx):
    t, r = relative_indices[idx]
    return f"{t}{r}"


# ============================================================
# PBC-AWARE GEOMETRY HELPERS
# ============================================================

def min_image_delta(i, j):
    d = coords[j] - coords[i]
    frac = d @ inv_lattice
    frac -= np.round(frac)
    return frac @ lattice


def min_image_distance(i, j):
    return np.linalg.norm(min_image_delta(i, j))


def bond_change(i, j, disp):
    v0 = min_image_delta(i, j)
    L0 = np.linalg.norm(v0)
    if L0 < 1e-9:
        return 0.0
    v1 = v0 + (disp[j] - disp[i])
    return np.linalg.norm(v1) - L0


def angle_change(a, c, b, disp):
    def _ang(v1, v2):
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-9 or n2 < 1e-9:
            return 0.0
        return np.arccos(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))

    v1i = min_image_delta(c, a)
    v2i = min_image_delta(c, b)
    v1d = v1i + (disp[a] - disp[c])
    v2d = v2i + (disp[b] - disp[c])
    return np.degrees(_ang(v1d, v2d) - _ang(v1i, v2i))


def rel_disp(atom_idx, center_idx, disp):
    return disp[atom_idx] - disp[center_idx]


def perp_disp(atom_idx, center_idx, disp):
    vec = min_image_delta(center_idx, atom_idx)
    L = np.linalg.norm(vec)
    if L < 1e-9:
        return 0.0
    rd = rel_disp(atom_idx, center_idx, disp)
    proj = np.dot(rd, vec) / L
    perp = rd - proj * (vec / L)
    return np.linalg.norm(perp)


def get_local_frame(center_idx, a_idx, b_idx):
    va = min_image_delta(center_idx, a_idx)
    vb = min_image_delta(center_idx, b_idx)
    na, nb = np.linalg.norm(va), np.linalg.norm(vb)
    if na < 1e-9 or nb < 1e-9:
        return None, None
    va_hat, vb_hat = va / na, vb / nb

    normal = np.cross(va_hat, vb_hat)
    normal_mag = np.linalg.norm(normal)
    if normal_mag < 1e-9:
        return None, None
    n_hat = normal / normal_mag

    bisector = va_hat + vb_hat
    b_mag = np.linalg.norm(bisector)
    if b_mag < 1e-9:
        return None, None
    b_hat = bisector / b_mag

    t_hat = np.cross(n_hat, b_hat)
    t_hat /= np.linalg.norm(t_hat)
    return n_hat, t_hat


def classify_pair_motion(disp_a, disp_b, n_hat, t_hat, mag_threshold):
    results = {}
    a_n, b_n = np.dot(disp_a, n_hat), np.dot(disp_b, n_hat)
    a_t, b_t = np.dot(disp_a, t_hat), np.dot(disp_b, t_hat)

    avg_n = np.mean([abs(a_n), abs(b_n)])
    avg_t = np.mean([abs(a_t), abs(b_t)])

    if avg_n > mag_threshold and a_n != 0 and b_n != 0 and np.sign(a_n) == np.sign(b_n):
        results['wagging'] = avg_n
    if avg_n > mag_threshold and a_n != 0 and b_n != 0 and np.sign(a_n) != np.sign(b_n):
        results['twisting'] = avg_n
    if avg_t > mag_threshold and a_t != 0 and b_t != 0 and np.sign(a_t) == np.sign(b_t):
        results['rocking'] = avg_t
    return results


# ============================================================
# NEIGHBOR TOPOLOGY
# ============================================================

def _compute_all_neighbors():
    cache = {}
    for i in range(num_atoms):
        neighbors = []
        for j in range(num_atoms):
            if j == i:
                continue
            cutoff = bond_cutoff(atom_labels[i], atom_labels[j])
            if cutoff is not None and min_image_distance(i, j) < cutoff:
                neighbors.append((j, atom_labels[j]))
        cache[i] = neighbors
    return cache


NEIGHBOR_CACHE = _compute_all_neighbors()


def get_neighbors(idx):
    return NEIGHBOR_CACHE[idx]


# ============================================================
# DISPLACEMENT / RAMAN DATA / EIGENVECTOR LOADING
# ============================================================

def load_mode_displacements(mode_num):
    path = os.path.join(DISPLACEMENTS_DIR, f"{mode_num}.txt")
    disp_list = []
    try:
        with open(path) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 6:
                    disp_list.append([float(parts[3]), float(parts[4]), float(parts[5])])
                elif len(parts) == 4:
                    disp_list.append([float(parts[1]), float(parts[2]), float(parts[3])])
    except Exception as e:
        return None, f"Mode {mode_num}: error loading {path}: {e}"

    if len(disp_list) == 0:
        return None, f"Mode {mode_num}: file is empty or unparseable ({path})"

    try:
        disp = np.array(disp_list).reshape(num_atoms, 3)
    except ValueError:
        return None, f"Mode {mode_num}: expected {num_atoms} atoms in {path}, got {len(disp_list)}"

    return disp, None


def load_raman_data(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Raman data file not found: {path}")
    modes_data = []
    warnings = []
    with open(path) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue
            p = line.split()
            if len(p) < 5:
                warnings.append(f"'{path}' line {line_num}: skipped (insufficient cols)")
                continue
            try:
                modes_data.append((int(p[0]), float(p[1]), float(p[4])))
            except ValueError:
                warnings.append(f"'{path}' line {line_num}: non-numeric data skipped")
    return modes_data, warnings


def load_eigenvectors(path):
    if not os.path.isfile(path):
        return {}
    eigenvectors = {}
    with open(path, errors="ignore") as f:
        for line in f:
            m = re.match(r'\s*(\d+)\s+f(?:/i)?\s*=\s*[\d.+\-Ee]+\s+THz', line)
            if not m or 'f/i' in line:
                continue
            mode_num = int(m.group(1))
            try:
                next(f)  # skip the "X Y Z dx dy dz" column header
            except StopIteration:
                continue
            dx = []
            for _ in range(num_atoms):
                try:
                    nxt = next(f)
                except StopIteration:
                    break
                parts = nxt.split()
                if len(parts) >= 6:
                    try:
                        dx.extend([float(parts[3]), float(parts[4]), float(parts[5])])
                    except ValueError:
                        break
            if len(dx) == num_atoms * 3:
                eigenvectors[mode_num] = np.array(dx).reshape(num_atoms, 3)

    return eigenvectors


def determine_symmetry(mode_num, eigenvectors):
    if mode_num not in eigenvectors:
        return "?", None
    b_vec = lattice[1]
    b_norm = np.linalg.norm(b_vec)
    if b_norm < 1e-12:
        return "?", None
    b_hat = b_vec / b_norm

    disp = eigenvectors[mode_num]
    b_proj = np.sum(disp * b_hat, axis=1)
    in_plane = disp - np.outer(b_proj, b_hat)

    in_plane_mag = np.linalg.norm(in_plane)
    out_of_plane_mag = np.linalg.norm(b_proj)
    total = np.sqrt(in_plane_mag ** 2 + out_of_plane_mag ** 2)
    if total < 1e-12:
        return "?", None

    ratio = out_of_plane_mag / total
    label = r"$A''$" if ratio > 0.5 else r"$A'$"
    return label, ratio


# ============================================================
# LOCAL MODE ANALYSIS
# ============================================================

def analyze_local_modes(disp, freq, thresholds):
    DISP_THR = thresholds['disp']
    BOND_THR = thresholds['stretch']
    ANGLE_THR = ANGLE_CHANGE_THRESHOLD

    results = []

    # 'unit' tags the magnitude scale: length motions are in Angstrom,
    # angle motions in degrees. Keeping them separate avoids ranking a
    # scissoring angle (degrees) against a stretch/wag amplitude (A).
    def add(desc, latex, mag, category, unit):
        results.append({'desc': desc, 'latex': latex, 'mag': mag,
                        'category': category, 'unit': unit})

    # --- Carbon-centered groups ---
    for ci in group_indices.get('C', []):
        nbrs = get_neighbors(ci)
        h_nbrs = [j for j, lbl in nbrs if lbl == 'H']
        cl_nbrs = [j for j, lbl in nbrs if lbl == 'Cl']
        n_nbrs = [j for j, lbl in nbrs if lbl == 'N']

        if len(h_nbrs) == 2:
            h1, h2 = h_nbrs
            d1, d2 = rel_disp(h1, ci, disp), rel_disp(h2, ci, disp)
            if np.linalg.norm(d1) > DISP_THR or np.linalg.norm(d2) > DISP_THR:
                b1, b2 = bond_change(ci, h1, disp), bond_change(ci, h2, disp)
                if freq > CH_STRETCH_MIN and abs(b1) > BOND_THR and abs(b2) > BOND_THR:
                    tag = "symmetric" if b1 * b2 > 0 else "asymmetric"
                    latex_tag = "s" if b1 * b2 > 0 else "as"
                    mag = np.mean([abs(b1), abs(b2)])
                    add(f"CH2 {tag} stretching ({fmt(ci)})",
                        rf"$\nu_{{\mathrm{{{latex_tag}}}}}$(CH$_2$)", mag, "CH_stretch", "length")

                ac = angle_change(h1, ci, h2, disp)
                if abs(ac) > ANGLE_THR:
                    add(f"CH2 scissoring ({fmt(h1)}-{fmt(ci)}-{fmt(h2)}, {ac:.2f} deg)",
                        r"$\beta$(CH$_2$)", abs(ac), "CH_bend", "angle")

                n_hat, t_hat = get_local_frame(ci, h1, h2)
                if n_hat is not None:
                    pm = classify_pair_motion(d1, d2, n_hat, t_hat, DISP_THR)
                    if 'wagging' in pm:
                        add(f"CH2 wagging ({fmt(h1)}, {fmt(h2)})", r"$\omega$(CH$_2$)", pm['wagging'], "CH_bend", "length")
                    if 'twisting' in pm:
                        add(f"CH2 twisting ({fmt(h1)}, {fmt(h2)})", r"$\tau$(CH$_2$)", pm['twisting'], "CH_bend", "length")
                    if 'rocking' in pm:
                        add(f"CH2 rocking ({fmt(h1)}, {fmt(h2)})", r"$\rho$(CH$_2$)", pm['rocking'], "CH_bend", "length")

        elif len(h_nbrs) == 3:
            h_disps = [rel_disp(h, ci, disp) for h in h_nbrs]
            if any(np.linalg.norm(d) > DISP_THR for d in h_disps):
                bcs = [bond_change(ci, h, disp) for h in h_nbrs]
                if freq > CH_STRETCH_MIN and all(abs(b) > BOND_THR for b in bcs):
                    symmetric = all(b > 0 for b in bcs) or all(b < 0 for b in bcs)
                    tag = "symmetric" if symmetric else "asymmetric"
                    latex_tag = "s" if symmetric else "as"
                    add(f"CH3 {tag} stretching ({fmt(ci)})",
                        rf"$\nu_{{\mathrm{{{latex_tag}}}}}$(CH$_3$)", np.mean(np.abs(bcs)), "CH_stretch", "length")

                umb = np.mean(h_disps, axis=0)
                if np.linalg.norm(umb) > DISP_THR:
                    add(f"CH3 umbrella ({fmt(ci)})", r"$\delta_{\mathrm{s}}$(CH$_3$)", np.linalg.norm(umb), "CH_bend", "length")

        for cl in cl_nbrs:
            bc = bond_change(ci, cl, disp)
            if abs(bc) > BOND_THR:
                dot = np.dot(disp[ci], disp[cl])
                tag = "symmetric" if dot > 0 else "asymmetric"
                latex_tag = "s" if dot > 0 else "as"
                add(f"C-Cl {tag} stretching ({fmt(ci)}-{fmt(cl)})",
                    rf"$\nu_{{\mathrm{{{latex_tag}}}}}$(C--Cl)", abs(bc), "framework_stretch", "length")

        for cl in cl_nbrs:
            for n in n_nbrs:
                ac = angle_change(cl, ci, n, disp)
                if abs(ac) > ANGLE_THR:
                    add(f"Cl-C-N bending ({fmt(cl)}-{fmt(ci)}-{fmt(n)})",
                        r"$\delta$(Cl--C--N)", abs(ac), "framework_bend", "angle")

    # --- Nitrogen centers ---
    for ni in group_indices.get('N', []):
        c_nbrs = [j for j, lbl in get_neighbors(ni) if lbl == 'C']
        if not c_nbrs:
            continue
        bcs = [bond_change(ni, c, disp) for c in c_nbrs]
        if any(abs(b) > BOND_THR for b in bcs):
            symmetric = all(b > 0 for b in bcs) or all(b < 0 for b in bcs)
            tag = "symmetric" if symmetric else "asymmetric"
            latex_tag = "s" if symmetric else "as"
            add(f"C-N {tag} stretching ({fmt(ni)})",
                rf"$\nu_{{\mathrm{{{latex_tag}}}}}$(C--N)", np.mean(np.abs(bcs)), "framework_stretch", "length")

    # --- Mn centers ---
    for zi in group_indices.get('Mn', []):
        cl_nbrs = [j for j, lbl in get_neighbors(zi) if lbl == 'Cl']
        if len(cl_nbrs) < 4:
            continue
        bcs = [bond_change(zi, c, disp) for c in cl_nbrs]
        if any(abs(b) > BOND_THR for b in bcs):
            symmetric = all(b > 0 for b in bcs) or all(b < 0 for b in bcs)
            tag = "symmetric" if symmetric else "asymmetric"
            latex_tag = "s" if symmetric else "as"
            add(f"Mn-Cl {tag} stretching ({fmt(zi)})",
                rf"$\nu_{{\mathrm{{{latex_tag}}}}}$(Mn--Cl)", np.mean(np.abs(bcs)), "framework_stretch", "length")

        max_ac = 0.0
        for i in range(len(cl_nbrs)):
            for j in range(i + 1, len(cl_nbrs)):
                max_ac = max(max_ac, abs(angle_change(cl_nbrs[i], zi, cl_nbrs[j], disp)))
        if max_ac > ANGLE_THR:
            add(f"MnCl6 deformation ({fmt(zi)})", r"$\delta$(MnCl$_6$)", max_ac, "framework_bend", "angle")

    # --- Cl-Cl interactions ---
    visited = set()
    for cl1 in group_indices.get('Cl', []):
        for cl2, lbl in get_neighbors(cl1):
            if lbl != 'Cl':
                continue
            key = tuple(sorted((cl1, cl2)))
            if key in visited:
                continue
            visited.add(key)
            cc = bond_change(cl1, cl2, disp)
            if abs(cc) > BOND_THR:
                add(f"Cl-Cl stretching ({fmt(cl1)}-{fmt(cl2)})", None, abs(cc), "framework_stretch", "length")

    length_res = sorted((r for r in results if r['unit'] == 'length'),
                        key=lambda x: x['mag'], reverse=True)
    angle_res = sorted((r for r in results if r['unit'] == 'angle'),
                       key=lambda x: x['mag'], reverse=True)
    return length_res + angle_res


# ============================================================
# MAIN EXECUTION ROUTINE
# ============================================================

def main():
    try:
        raman_modes, warnings = load_raman_data(RAMAN_FILE)
        for w in warnings:
            print(f"WARNING: {w}", file=sys.stderr)
    except Exception as e:
        print(f"ERROR reading Raman data: {e}", file=sys.stderr)
        sys.exit(1)

    eigenvectors = load_eigenvectors(OUTCAR_PHON_FILE)

    report_lines = []
    csv_rows = []
    latex_rows = []

    for mode_num, freq, activity in raman_modes:
        if not (FREQ_MIN <= freq <= FREQ_MAX):
            continue

        disp, err = load_mode_displacements(mode_num)
        if err:
            print(f"WARNING: {err}", file=sys.stderr)
            continue

        max_disp = np.max(np.linalg.norm(disp, axis=1))
        thresholds = {
            'disp': max_disp * DISP_ADAPTIVE_FRACTION,
            'stretch': max(BOND_ABS_THRESHOLD, max_disp * DISP_ADAPTIVE_FRACTION)
        }

        mode_results = analyze_local_modes(disp, freq, thresholds)
        symmetry_label, _ = determine_symmetry(mode_num, eigenvectors)

        unit_suffix = lambda u: "\u00c5" if u == "length" else "\u00b0"

        # Build Full Report Output
        report_lines.append(f"Mode {mode_num:3d} | Freq: {freq:7.2f} cm^-1 | Activity: {activity:8.3f}")
        n_logged = 0
        for res in mode_results:
            report_lines.append(
                f"  - [{res['category']:<16}] {res['desc']} "
                f"(mag={res['mag']:.4f} {unit_suffix(res['unit'])})"
            )
            n_logged += 1
            if n_logged >= REPORT_TOP_N:
                break
        report_lines.append("")

        # Build CSV Rows
        top_descs = [r['desc'] for r in mode_results[:CSV_TOP_N]]
        csv_rows.append([mode_num, f"{freq:.2f}", f"{activity:.3f}", symmetry_label] + top_descs)

        # Build LaTeX Rows: length motions compete via the 0.30 relative
        # filter; angle motions fill remaining slots by |angle change| so
        # degrees never crowd out Angstrom-scale assignments.
        latex_labels = []
        length_res = [r for r in mode_results if r['unit'] == 'length']
        angle_res = [r for r in mode_results if r['unit'] == 'angle']
        if length_res:
            top_mag = length_res[0]['mag']
            for res in length_res:
                if res['latex'] and (res['mag'] >= top_mag * RELATIVE_ASSIGNMENT_THRESHOLD):
                    if res['latex'] not in latex_labels:
                        latex_labels.append(res['latex'])
        for res in angle_res:
            if res['latex'] and res['latex'] not in latex_labels:
                latex_labels.append(res['latex'])
            if len(latex_labels) >= MAX_LATEX_ASSIGNMENTS:
                break
        assignment_str = ", ".join(latex_labels) if latex_labels else "Unassigned"
        latex_rows.append(f"{mode_num} & {freq:.1f} & {symmetry_label} & {assignment_str} \\\\")

    # Output Files Writing
    with open(OUTPUT_REPORT_FILE, "w") as f:
        f.write("\n".join(report_lines))

    with open(OUTPUT_SUMMARY_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Mode", "Frequency_cm-1", "Raman_Activity", "Symmetry", "Assign_1", "Assign_2", "Assign_3", "Assign_4"])
        writer.writerows(csv_rows)

    with open(OUTPUT_LATEX_FILE, "w") as f:
        f.write("% Raman Table Output\n")
        f.write("\\begin{tabular}{cccc}\n\\hline\nMode & Freq (cm$^{-1}$) & Sym & Assignment \\\\\n\\hline\n")
        f.write("\n".join(latex_rows))
        f.write("\n\\hline\n\\end{tabular}\n")

    print(f"Analysis Complete.")
    print(f"  Report -> {OUTPUT_REPORT_FILE}")
    print(f"  CSV    -> {OUTPUT_SUMMARY_CSV}")
    print(f"  LaTeX  -> {OUTPUT_LATEX_FILE}")


if __name__ == "__main__":
    main()