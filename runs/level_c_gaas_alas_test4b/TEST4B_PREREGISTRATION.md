# TEST 4B PREREGISTRATION — FC1 GaAs/AlAs [001] interface DFT and D_AB

Written before any FC1 force or IFC evaluation. All thresholds below are frozen.

## 1. Frozen FC1 definition (Test 4, unchanged)
`{"name": "FC1_ideal_abrupt", "description": "chemically abrupt coherent interface; in-plane lattice = GaAs (substrate-constrained), AlAs tetragonally strained", "a_par_A": 5.658372605963021, "a_perp_GaAs_A": 5.658372605963021, "a_perp_AlAs_A": 5.689594201107272, "intermixing": null, "registry": "As-terminated shared anion plane; Ga plane | As plane | Al plane", "supercell": "[001] SL cell, 1x1 conventional in-plane, >=6 ML each side; periodic; 2 interfaces per cell (one is the mirror)", "eps_par_AlAs": -0.0028171735699556726, "status": "DEFINED (structure/registry/strain frozen); interface IFCs NOT COMPUTED - requires DFT (no electronic-structure code available on this platform)", "ifc_status": "NOT COMPUTATIONALLY FEASIBLE ON THIS PLATFORM (no DFT code); D_AB absent"}`
Test 4 manifest sha256 `98d8590eebc7362386e3a6617941ded33217ee715dec55ba5092985b95683013`; FC1 definition.json sha256 `d45484a2f827b71e0112f60564849f34720a9a8e1545a075faac0493c1226c01`.

## 2. Supercell
nA = nB = 8 monolayer pairs (4a per side ≥ 2 × IFC range 2a), 64 atoms, 1×1 conventional in-plane (frozen FC1), c = 45.3919 Å, 2 symmetry-equivalent As-centred interfaces (−4 axis).
Structure hash (cell, species, fractional coordinates): `dd413261893b9ae477c49c7d9ea0d20cc831f1e494d3f249a80057269b1d04c0`. Initial files: `structures/FC1_initial_POSCAR`, `structures/FC1_initial_structure.json`.
Technical note (not a definition change): Test 4 wrote "≥ 6 ML each side"; 8 pairs per side is chosen so that ≥ 2a of bulk-like material separates the two periodic interfaces from the IFC range of each (bulk fc2 |Φ| < 0.04 eV/Å² beyond 9.5 Å).

## 3. DFT environment
Quantum ESPRESSO pw.x 7.5 (conda-forge build, Docker image rcm-qe:test4b); pseudopotentials: Ga: Ga.pbesol-dn-kjpaw_psl.1.0.0.UPF (sha256 60b44b95c5547bb7…), As: As.pbesol-n-kjpaw_psl.1.0.0.UPF (sha256 0bbf3763c2095dbd…), Al: Al.pbesol-n-kjpaw_psl.1.0.0.UPF (sha256 0a36cbca2d703e17…)

## 4. Settings and deviations from the bulk methodology
- code: VASP -> QE (no VASP licence/executable on this machine); same plane-wave PAW/PBEsol methodology
- PAW datasets: VASP PAW_PBE Ga_d/As_d/Al -> PSlibrary 1.0.0 PBEsol Ga (3d in valence), As (3d frozen; no As_d PBEsol PAW in PSlibrary), Al
- cutoffs: ENCUT 520 eV -> ecutwfc 60 Ry (816 eV) / ecutrho 480 Ry, the PSlibrary-recommended minimum for Ga-d being 60 Ry
- occupations: ISMEAR=0, SIGMA=0.01 eV -> fixed (insulating, equivalent)
- k-points: 4x4x1 shifted Monkhorst-Pack for the 1x1x(16 pairs) SL cell, matching the archived force-cell density (shifted 1x1x1 on 4a)
- settings: {"ecutwfc_Ry": 60.0, "ecutrho_Ry": 480.0, "conv_thr_Ry": 1e-10, "occupations": "fixed", "mixing_beta": 0.4, "forc_conv_thr_eV_A": 0.0001, "etot_conv_thr_Ry": 1e-07, "press_conv_thr_kbar": 0.1, "displacement_A": 0.03}; k-points (4, 4, 1) shift (1, 1, 0); relaxation vc-relax with cell_dofree='z'.
- Bulk cross-validation (gate BX, must pass before FC1 forces are used): GaAs and AlAs: vc-relax (8x8x8 shifted MP, primitive) -> a; FD (0.03 A) in 2x2x2 conventional cells (2x2x2 shifted MP) -> Gamma/X frequencies (no NAC) and k_par=0 layer blocks vs Test 4; thresholds {"a_rel": 0.005, "gamma_TO_noNAC_rel": 0.03, "X_acoustic_rel": 0.03, "layer_block_H0_H1_rel_frobenius": 0.05, "note": "QE/PSlibrary vs Test 4 VASP-derived bulk (same PBEsol functional): lattice constant from vc-relax; frequencies from a 2x2x2 conventional FD supercell without NAC; [001] layer blocks at k_par=0"}

## 5. Gates T4B-1 … T4B-10 (frozen thresholds)
- **T4B_1_structure**: {"rule": "FC1 cell built from the frozen Test 4 definition; structure hash recorded before any calculation"}
- **T4B_2_relaxation**: {"F_max_eV_A": 0.0001, "sigma_zz_kbar": 0.1, "bond_length_rel_to_bulk": 0.03, "max_displacement_A": 0.15, "c_change_rel": 0.01, "rule": "vc-relax with cell_dofree=z (in-plane fixed = coherent to GaAs); positions free; FC1 identity preserved: no interplanar swap, no bond broken/created, both interfaces remain As-centred Ga|As|Al"}
- **T4B_3_ifc**: {"displacement_A": 0.03, "asr_raw_eV_A2": 0.01, "asr_corrected_eV_A2": 1e-08, "conv_check_forces_rel": 0.02, "rule": "phonopy symmetry-reduced finite displacements on the relaxed cell; pw.x scf with identical settings (nosym); fc2 by phonopy; symmetrization = translational ASR + permutation; raw and corrected residuals both reported"}
- **T4B_4_DAB**: {"rule": "D_AB = mass-weighted fc2 block coupling GaAs-side atoms (pairs < nA) to AlAs-side atoms (pairs >= nA) at k_par = 0 (1x1 in-plane cell = in-plane image sum), units eV/A^2/amu; atom indices, planes, Cartesian ordering and range recorded; block symmetry D_BA = D_AB^T checked"}
- **T4B_5_bulk_recovery**: {"central_block_rel_frobenius": 0.05, "decay_monotone": true, "rule": "per-plane on-site (3x3 summed over the plane) and nearest-plane coupling blocks of the SL compared with the Test 4 bulk H0/H1 blocks at k_par=0; deviation at the central planes (>= 2a from any interface) <= 5 % (relative Frobenius)"}
- **T4B_6_stability**: {"negative_freq_tolerance_THz": -0.05, "rule": "SL dynamical matrix at k_par=0 for q_z on 8 points; no eigenfrequency below -0.05 THz except the 3 translational zero modes; localized modes reported with plane participation, not rejected"}
- **T4B_7_agf_integration**: {"rule": "device = 4 GaAs pairs + 4 AlAs pairs (96 DOF) cut from the SL at the slab centres; leads = Test 4 PL blocks (hash-checked); device-lead coupling = bulk H1/H2 (justified by T4B-5); Test 4 agf3d code unchanged"}
- **T4B_8_agf_conservation**: {"positivity": -1e-09, "channel_bound_tol": 0.005, "unitarity": 0.005, "reciprocity": 1e-05, "eta_rel": [0.0001, 1e-05, 1e-06], "sgf_tol": [1e-12, 1e-14], "edge_window_THz": 0.02, "rule": "same criteria as Test 4 smoke test, evaluated on the FC1 interface device; eta and SGF-tolerance convergence; frequency-grid independence (50 vs 200 points)"}
- **T4B_9_continuum**: {"band_THz": [0.05, 0.4], "abs_tol_total": 0.02, "abs_tol_L": 0.02, "rule": "T_FC1(omega) in the long-wavelength band vs matched continuum: total = T_L^c + 2 T_T^c with T^c = 4 Z_A Z_B/(Z_A+Z_B)^2 per polarization (Test 4 parameters); L channel via z-polarized projection; no tuning"}
- **T4B_10_provenance**: {"rule": "every D_AB element traces to: pw.x outputs (hashed) -> phonopy displacement dataset (hashed) -> fc2 (hashed) -> block extraction script (hashed)"}

## 6. IFC convergence assessments (Sec. 10)
- thickness: bulk recovery (T4B-5) across the slab; - in-plane size: 1×1 = frozen FC1; the k∥=0 blocks are exact for this cell, k∥≠0 blocks are NOT produced (limitation stated, no embodiment change); - displacement amplitude: ± pair on the interface As (force antisymmetry ≤ 2 %); - electronic convergence: conv_thr 1e-10 vs 1e-8 on one displaced cell (ΔF ≤ 1e-4 eV/Å); - cutoff: bulk FD forces at 50 vs 60 Ry (≤ 2 %); - symmetry: phonopy/spglib symprec 1e-5, displaced cells run with nosym.

## 7. Prospective Series 3B mapping rule (Sec. 19)
- unit: one [001] monolayer pair u_i = a_perp(m_i)/2 in the FC1 strain state (GaAs 0.28292 nm, AlAs 0.28448 nm)
- allocation: n_i^0 = max(1, round(t_i/u_i)); DeltaL = sum n_i^0 u_i - L_target; m = round(|DeltaL|/u_bar); move the m layers with the largest residual (n_i^0 u_i - t_i) of the correcting sign by one unit each (each layer at most once, ties -> lowest index)
- guarantees: |Delta t_i| <= u_i for every layer (rounding <= u_i/2, at most one correcting unit); |Delta L| <= u_bar (one unit)
- acceptance: {"per_layer_ratio_max": 1.0, "rms_delta_t_nm": 0.1, "abs_delta_L_nm": 0.15, "sub_unit_layers": "count reported (t_i < u_i mapped to n_i = 1), not gated"}
- derivation: per-layer bound = discrete realizability (one repeat unit); RMS bound = uniform rounding RMS u/sqrt(12) = 0.082 nm plus correction moves (<= 0.10 nm); total-thickness bound = half a unit rounding of the global count plus material-unit mismatch (<= 0.15 nm)

## 8. Feasibility rule
if the projected FC1 workflow (relaxation + all displacements) exceeds 72 h of wall time on this machine, execute what fits and classify NOT READY — FC1 DFT REQUIRED; never reduce k-points/cutoffs below the preregistered values to fit

## 9. Forbidden
No Series 3B, no sensing metrics, no geometry optimisation, no fitted/averaged/impedance-derived interface constants, no parameter changes after inspecting phonon/transport results.