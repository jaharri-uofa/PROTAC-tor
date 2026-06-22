# PROTACtor Scoring & Analysis Documentation

This document explains all scores and metrics calculated during PROTACtor analysis, organized by pipeline stage.

---

## Table of Contents

1. [Protein–Protein Docking](#protein--protein-docking)
2. [Linker Generation](#linker-generation)
3. [PROTAC Docking](#protac-docking)
4. [Molecular Dynamics (MD) Analysis](#molecular-dynamics-md-analysis)
5. [Advanced Metrics](#advanced-metrics)

---

## Protein–Protein Docking

### Compatibility Score

**Description:**  
Assesses how well the ligase (E3) and POI (target protein) interact in the docked complex. Calculated by averaging the docking scores from all generated complexes (those with valid scores) and normalizing by a reference value.

**Formula:**
```
Compatibility Score = (Sum of all complex scores) / (Number of complexes) / 1990
```

**Interpretation:**
- Higher scores indicate better ligase-POI compatibility
- Scores closer to 1.0 suggest good protein-protein compatibility
- Lower scores may indicate poor geometric fit or unfavorable interactions

**Source Files:**
- `complex.*.pdb` (generated during protein-protein docking)

**Units:** Normalized dimensionless score

---

### Ligand Distance Range

**Description:**  
Reports the minimum and maximum distances (in Ångströms) between ligand atoms and the docking reference point throughout the docking stage.

**Calculation:**
- Parses `lig_distances.txt` for all reported distances
- Extracts numeric values from lines containing a colon
- Returns min and max values

**Interpretation:**
- **Minimum distance:** Closest approach of any ligand atom
- **Maximum distance:** Farthest ligand atom position
- Useful for assessing ligand positioning consistency

**Source Files:**
- `lig_distances.txt`

**Units:** Ångströms (Å)

---

## Linker Generation

### Warhead SMILES

**Description:**  
Chemical structure representation of the two warheads (bioactive molecules) connected by the generated linkers. SMILES notation allows unambiguous representation of molecular structure.

**Interpretation:**
- **Warhead 1 (POI-targeting):** The primary ligand targeting your protein of interest
- **Warhead 2 (E3-targeting):** The ligand targeting the E3 ligase
- Both warheads are held together by optimized linker chains

**Source Files:**
- `smiles.smi`

---

### Total Linkers Sampled

**Description:**  
The total number of unique linker structures generated and evaluated during linker generation (Stage 1 of LinkInvent).

**Interpretation:**
- Higher numbers indicate more extensive sampling
- Provides a measure of chemical space exploration
- Useful for assessing optimization convergence

**Source Files:**
- `linkinvent_stage_1.csv`

---

### Best Linker Score

**Description:**  
The highest (best) linkage score achieved during linker generation, indicating the quality of the top-ranked optimized linker.

**Interpretation:**
- Higher is better
- Represents the model's confidence in the generated linker
- Scores typically range from 0 to 1 or higher depending on the scoring model
- Guides selection of top candidates for docking

**Source Files:**
- `linkinvent_stage_1.csv` (Score column, maximum value)

---

## PROTAC Docking

### Warhead Binding Affinity

**Description:**  
Average binding affinity of the warheads to their respective targets, calculated from the top-ranked PROTAC complexes.

**Calculation:**
- Parses `top5_complexes.csv` for affinity values
- Extracts numeric affinity values (typically in kcal/mol)
- Computes mean and standard deviation

**Interpretation:**
- **More negative values = stronger binding**
- Typical range: -5 to -15 kcal/mol
- Represents independent binding strength of the warheads
  - useful for comparison to PROTAC binding affinity
- Standard deviation indicates consistency across top candidates

**Units:** kcal/mol

**Source Files:**
- `top5_complexes.csv`

---

### Highest PROTAC Affinity

**Description:**  
The single best (most negative/strongest) binding affinity score from all PROTAC ternary complexes evaluated during docking.

**Interpretation:**
- Represents the "best-case" binding scenario
- More negative = stronger binding
- Used to identify the optimal PROTAC candidate
- Lower (more negative) values indicate tighter complex stability

**Units:** kcal/mol (typically negative)

**Source Files:**
- `top5_complexes.csv`

---

### Top PROTAC Complexes Table

**Description:**  
Ranked list of the top 5 PROTAC ternary complexes by binding affinity, showing their SMILES structure and calculated binding free energy.

**Columns:**
- **Rank:** Position in sorted order (1 = best/lowest affinity)
- **SMILES:** Chemical structure notation of the complete PROTAC molecule
- **Affinity:** Predicted binding free energy (more negative = better)

**Interpretation:**
- Guides experimental validation prioritization
- Shows chemical diversity of top candidates
- Affinity spread indicates robustness of design

**Source Files:**
- `top5_complexes.csv`

---

## Molecular Dynamics (MD) Analysis

MD analysis evaluates complex stability and dynamics through trajectory analysis of production MD runs. Two types of systems are analyzed:

1. **PROTAC Ternary Complex:** E3 ligase + POI + PROTAC linker
2. **Control (Warhead-only):** E3 ligase + POI + warhead (no linker)

### Root Mean Square Deviation (RMSD) — Backbone

**Description:**  
Measures how much the protein backbone (Cα, C, N atoms only) deviates from its initial structure over the course of the MD simulation.

**Calculation:**
- Uses `cpptraj` to compute RMSD vs. first frame
- Measures distance between corresponding backbone atoms over time
- Includes water and ions stripped out for clarity

**Formula:**
```
RMSD(t) = sqrt( (1/N) * Σ(distance_i(t))² )
```
where N = number of atoms and distance_i = distance of atom i from initial position

**Interpretation:**
- **Equilibration phase:** Initial rise in RMSD as system relaxes
- **Plateau/stable:** System has reached stable conformation
- **Higher values:** More structural fluctuation
- **Ideal:** RMSD stabilizes and remains relatively constant
- **> 2-3 Å:** May indicate instability or improper solvation

**Statistics Reported:**
- **Mean RMSD:** Average over entire trajectory (indicates overall stability)
- **Std Dev:** Standard deviation of fluctuations around mean

**Units:** Ångströms (Å)

**Output Files:**
- `*_backbone_rmsd.png` (plot)
- `*_backbone_rmsd.dat` (raw data)

---

### Root Mean Square Deviation (RMSD) — Ligand

**Description:**  
Measures positional deviation of the ligand (PROTAC or warhead) from its initial docked position, excluding hydrogen atoms.

**Calculation:**
- Heavy atoms only (non-hydrogen)
- No fit mode (absolute positions relative to frame 1)
- Computed using PROTAC residue name identified from ligand_resname.txt

**Interpretation:**
- **Low RMSD:** Ligand remains in docked pose (good binding stability)
- **High RMSD:** Ligand displaces or dissociates
- **Spike then stable:** Initial repositioning followed by binding
- **Continuously rising:** Indicates ligand unbinding

**Thresholds:**
- **< 1.5 Å:** Excellent stability
- **1.5–3.0 Å:** Good, minor displacement
- **> 3.0 Å:** Significant movement, potential instability

**Units:** Ångströms (Å)

**Output Files:**
- `*_ligand_rmsd.png` (plot)
- `*_ligand_rmsd.dat` (raw data)

---

### Root Mean Square Fluctuation (RMSF) — Per-Residue

**Description:**  
Measures the average displacement of each individual residue's backbone (Cα, C, N) from its mean position over the entire simulation.

**Calculation:**
- For each residue, computes average backbone atom position
- Then calculates deviation of each frame from that mean
- Provides per-residue flexibility profile

**Formula:**
```
RMSF_residue_i = sqrt( (1/T) * Σ_t(distance_i,t)² )
```
where T = number of trajectory frames

**Interpretation:**
- **Low RMSF (< 1 Å):** Rigid, well-ordered region
- **Moderate RMSF (1–2 Å):** Typical for structured regions
- **High RMSF (> 3 Å):** Flexible loop or dynamic region
- **Peaks in RMSF:** Identify functionally important mobile regions

**Use Cases:**
- Identify conformational hot-spots
- Assess linker-induced dynamics
- Compare PROTAC vs. control flexibility

**Statistics Reported:**
- **Max RMSF residue:** Residue number with highest flexibility
- **Max RMSF value:** Magnitude of maximum fluctuation

**Units:** Ångströms (Å)

**Output Files:**
- `*_rmsf.png` (per-residue bar/line plot)
- `*_rmsf_byres.dat` (raw data)

---

## Advanced Metrics

### MM/GBSA: Binding Free Energy (ΔG)

**Description:**  
Molecular Mechanics - Generalized Born Surface Area (MM/GBSA) calculation estimates the absolute binding free energy of the ternary complex using molecular force fields.

**Components:**
```
ΔG_binding = ΔG_gas + ΔG_solvation - T*ΔS

where:
  ΔG_gas          = ΔE_internal + ΔE_electrostatic + ΔE_VdW
  ΔG_solvation    = ΔG_polar + ΔG_nonpolar
  T*ΔS           = entropic contribution (often approximated)
```

**Interpretation:**
- **Negative ΔG:** Favorable/spontaneous binding (more negative = stronger)
- **Positive ΔG:** Unfavorable binding
- **Typical range:** -40 to -80 kcal/mol for specific interactions
  - **NOTE:** the MMGBSA values are inflated due to ignoring entropy contributions, use as qualtitaive data ONLY
- More reliable than pure docking scores due to MD sampling

**Output Nomenclature:**
- `out_ligand.dat` or `out_ligand1.dat`: Primary PROTAC/warhead binding energy
- `out_ligand2.dat` (if multi-ligand system): Secondary component

**Limitations:**
- Entropy term (T*ΔS) often approximated or ignored
- Assumes MM force field accuracy
- Sensitive to solvation model choice
- Best used for relative comparisons (same system, different conditions)

**Units:** kcal/mol

**Output Files:**
- `*_out_ligand*.dat` (parsed statistics)

---

### ΔΔG: Linker Contribution to Binding

**Description:**  
Differential binding free energy quantifying the energetic contribution of the linker to ternary complex stabilization, computed by comparing PROTAC and control (warhead-only) systems.

**Formula:**
```
ΔΔG_linker = ΔG(PROTAC ternary) − ΔG(control warheads only)
           = ΔG_PROTAC − ΔG_control
```

**Interpretation:**
- **Negative ΔΔG:** Linker stabilizes the ternary complex ✓ (GOOD)
  - The linker energy contribution is favorable
  - Ternary complex is more stable than warheads alone
  
- **Positive ΔΔG:** Linker destabilizes the ternary complex ✗ (BAD)
  - The linker imposes energetic strain
  - Warheads bind better independently

- **|ΔΔG| magnitude:**
  - Larger magnitude (more negative) = stronger contribution
  - Typical optimal range: -2 to -8 kcal/mol

**Comparison Interpretation:**
| ΔΔG Value | Linker Assessment |
|-----------|------------------|
| < -5 kcal/mol | Excellent stabilization |
| -2 to -5 kcal/mol | Good stabilization |
| -1 to 0 kcal/mol | Minimal/neutral effect |
| 0 to +2 kcal/mol | Minor destabilization |
| > +2 kcal/mol | Significant destabilization |

**Uncertainty:**
- Reported as ΔΔG ± std (combined error from both calculations)
- Propagated from PROTAC and control MM/GBSA uncertainties

**Conditions for Calculation:**
- Requires completed control MD run
- Requires at least one completed PROTAC MD run
- Both must have valid MM/GBSA calculations

**Source Files:**
- Control: `ternary_complex_control/md/out_ligand.dat`
- PROTAC: `ternary_complex*/md/out_ligand*.dat`

**Units:** kcal/mol

---

## Summary of Key Metrics

| Stage | Metric | Unit | Ideal Range | Interpretation |
|-------|--------|------|-------------|-----------------|
| **PP Docking** | Compatibility | Norm. | 0.7–1.0 | Higher = better E3–POI fit |
| **Linker Gen** | Best Score | Score | > 0.7 | Higher = better linker quality |
| **PROTAC Dock** | Affinity | kcal/mol | < -10 | More negative = stronger |
| **MD (BB RMSD)** | Mean RMSD | Å | 1–3 | Stable structure |
| **MD (Lig RMSD)** | Mean RMSD | Å | < 2 | Ligand stays bound |
| **MD (RMSF)** | Max per-residue | Å | < 3 | Most regions stable |
| **MM/GBSA** | ΔG | kcal/mol | -5 to -15 | More negative = favorable |
| **ΔΔG** | Linker energy | kcal/mol | < -2 | Negative = stabilizing linker |

---

## References

- **RMSD/RMSF:** Computed using AMBER `cpptraj` tool
- **MM/GBSA:** MMPBSA.py (from AMBER MD package)
- **Docking:** ZDOCK for rigid-body docking; SMILES via RDKit
- **Linker Optimization:** LinkInvent (machine learning model)

---

## Example Report Output

```
===================================================================
PROTACtor Analysis Report
===================================================================

## Protein–Protein Docking

  Compatibility score          : 0.8234
  Ligand distance range        : 2.45 – 18.92 Å

## Linker Generation

  Warhead SMILES (1)           : c1ccc(cc1)C(=O)O
  Warhead SMILES (2)           : Cc1ccc(cc1)S(=O)(=O)N
  Total linkers sampled        : 1250
  Best linker score            : 0.9123

## PROTAC Docking

  Warhead binding affinity     : -11.34 ± 1.20 kcal/mol
  Highest PROTAC affinity      : -13.45

  Top PROTAC complexes:
    Rank  smiles                                    affinity
    1     CC(=O)Nc1ccc(cc1)C(c2ccc...)          -13.45
    2     CC(=O)Nc1ccc(cc1)C(c2ccc...)          -13.12
    3     CC(=O)Nc1ccc(cc1)C(c2ccc...)          -12.87
    ...

## MD Analysis

### PROTAC Ternary Complexes

  Complex          : corrected_CDK5_ternary_complex1
  Backbone RMSD    : 2.15 ± 0.45 Å
  Ligand RMSD      : 1.82 ± 0.38 Å
  Max RMSF residue : residue 156 (2.89 Å)
  out_ligand ΔG   : -8.67 ± 1.23 kcal/mol

  Average PROTAC ΔG (all complexes) : -8.45 kcal/mol

### Control (Warheads Only)

  Backbone RMSD    : 1.98 ± 0.42 Å
  out_ligand ΔG   : -6.23 ± 1.10 kcal/mol

### Linker Contribution: ΔΔG = ΔG(PROTAC) − ΔG(control)

  Interpretation: negative ΔΔG = linker stabilises the ternary
  complex relative to the warhead-only control.
  
  PROTAC ΔG                         : -8.45 ± 0.95 kcal/mol
  Control ΔG                        : -6.23 ± 1.10 kcal/mol
  ΔΔG (linker contribution)         : -2.22 ± 1.45 kcal/mol
```

---

## Contact & Contributions

For questions about scoring methodology or to suggest improvements to the analysis pipeline, please refer to the main PROTACtor repository or contact the development team.
