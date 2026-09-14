"""Quantum ESPRESSO (pw.x) driver: input writer, output parser, containerised execution (Docker image rcm-qe:test4b).
All inputs/outputs are written under runs/level_c_gaas_alas_test4b/dft/<calc>/ and kept for reproduction."""
from __future__ import annotations
import os, re, subprocess, time, hashlib, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DFT = ROOT / 'runs/level_c_gaas_alas_test4b/dft'
IMAGE = 'rcm-qe:test4b'
PSEUDO = {'Ga': 'Ga.pbesol-dn-kjpaw_psl.1.0.0.UPF', 'As': 'As.pbesol-n-kjpaw_psl.1.0.0.UPF', 'Al': 'Al.pbesol-n-kjpaw_psl.1.0.0.UPF'}
MASSES = dict(Ga=69.723, As=74.9216, Al=26.9815385)
RY_EV = 13.605693122994; BOHR_A = 0.529177210903
RYBOHR_EVA = RY_EV / BOHR_A                                   # Ry/bohr -> eV/Å

# Preregistered settings (TEST4B_PREREGISTRATION.md Sec. 4) — do not change after FC1 forces are evaluated
SETTINGS = dict(ecutwfc_Ry=60.0, ecutrho_Ry=480.0, conv_thr_Ry=1e-10, occupations='fixed', mixing_beta=0.4,
                forc_conv_thr_eV_A=1e-4, etot_conv_thr_Ry=1e-7, press_conv_thr_kbar=0.1, displacement_A=0.03)


def write_pwin(path, calc, cell_A, symbols, frac, kpts=(4, 4, 1), kshift=(1, 1, 0), prefix='pw', cell_dofree=None, nosym=False, nbnd=None, extra_system=None, settings=SETTINGS, disk_io='medium', diagonalization='david'):
    """disk_io='medium' keeps one k-point's wavefunctions in memory (others on container-local scratch); diagonalization
    selects the iterative eigensolver — both are memory/algorithm choices with no effect on converged results."""
    species = []; [species.append(s) for s in symbols if s not in species]
    L = [f"&CONTROL\n  calculation='{calc}'\n  prefix='{prefix}'\n  pseudo_dir='/work/pseudo'\n  outdir='/tmp/qe_{prefix}'\n  tprnfor=.true.\n  tstress=.true.\n  disk_io='{disk_io}'\n  verbosity='low'\n"
         f"  forc_conv_thr={settings['forc_conv_thr_eV_A'] / RYBOHR_EVA:.3e}\n  etot_conv_thr={settings['etot_conv_thr_Ry']:.1e}\n/\n"
         f"&SYSTEM\n  ibrav=0\n  nat={len(symbols)}\n  ntyp={len(species)}\n  ecutwfc={settings['ecutwfc_Ry']}\n  ecutrho={settings['ecutrho_Ry']}\n  occupations='{settings['occupations']}'\n"
         + ("  nosym=.true.\n  noinv=.true.\n" if nosym else '') + (f"  nbnd={nbnd}\n" if nbnd else '') + (extra_system or '') + "/\n"
         f"&ELECTRONS\n  conv_thr={settings['conv_thr_Ry']:.1e}\n  mixing_beta={settings['mixing_beta']}\n  electron_maxstep=200\n  diagonalization='{diagonalization}'\n/\n"]
    if calc in ('relax', 'vc-relax'):
        L.append("&IONS\n  ion_dynamics='bfgs'\n/\n")
    if calc == 'vc-relax':
        L.append(f"&CELL\n  cell_dynamics='bfgs'\n  press_conv_thr={settings['press_conv_thr_kbar']}\n" + (f"  cell_dofree='{cell_dofree}'\n" if cell_dofree else '') + "/\n")
    L.append('ATOMIC_SPECIES\n' + ''.join(f'  {s} {MASSES[s]} {PSEUDO[s]}\n' for s in species))
    L.append('CELL_PARAMETERS angstrom\n' + ''.join('  ' + ' '.join(f'{v:.10f}' for v in row) + '\n' for row in np.asarray(cell_A)))
    L.append('ATOMIC_POSITIONS crystal\n' + ''.join(f'  {s} ' + ' '.join(f'{v:.10f}' for v in f) + '\n' for s, f in zip(symbols, np.asarray(frac))))
    L.append('K_POINTS automatic\n  ' + ' '.join(map(str, kpts)) + ' ' + ' '.join(map(str, kshift)) + '\n')
    Path(path).write_text(''.join(L)); return Path(path)


def run_pw(calc_dir: Path, np_=6, npool=1, timeout_s=48 * 3600, infile='pw.in', outfile='pw.out'):
    """Run pw.x inside the container with the calc directory bind-mounted. Returns wall time (s)."""
    calc_dir = Path(calc_dir); rel = calc_dir.relative_to(DFT).as_posix()
    cmd = ['docker', 'run', '--rm', '-v', f'{DFT.as_posix()}:/work', IMAGE, 'bash', '-lc',
           f'cd /work/{rel} && mpirun --oversubscribe --bind-to none -np {np_} pw.x -nk {npool} -in {infile} > {outfile} 2> pw.err; rm -rf /tmp/qe_*']
    t0 = time.time(); subprocess.run(cmd, check=False, timeout=timeout_s)
    err = (calc_dir / 'pw.err').read_text(errors='ignore') if (calc_dir / 'pw.err').exists() else ''
    if 'Killed' in err and np_ > 6:                                   # OOM under the 8 GB Docker cap: retry with fewer ranks (results unaffected)
        (calc_dir / 'pw.err.oom').write_text(err); return run_pw(calc_dir, np_=6, npool=1, timeout_s=timeout_s, infile=infile, outfile=outfile)
    return time.time() - t0


def parse_pwout(path) -> dict:
    txt = Path(path).read_text(errors='ignore')
    E = [float(x) for x in re.findall(r'^!\s+total energy\s+=\s+(-?\d+\.\d+) Ry', txt, flags=re.M)]
    forces = []
    for blk in re.findall(r'Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n\n((?:\s+atom\s+\d+ type\s+\d+\s+force =\s+\S+\s+\S+\s+\S+\n)+)', txt):
        forces.append(np.array([[float(v) for v in l.split('=')[1].split()] for l in blk.strip().splitlines()]) * RYBOHR_EVA)
    stress = []
    for blk in re.findall(r'total\s+stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*(-?\d+\.\d+)\n((?:.*\n){3})', txt):
        stress.append(np.array([[float(v) for v in l.split()[3:6]] for l in blk[1].strip().splitlines()]))     # kbar
    fin = re.search(r'Begin final coordinates(.*?)End final coordinates', txt, flags=re.S)
    cell = frac = syms = None
    if fin:
        s = fin.group(1); cm = re.search(r'CELL_PARAMETERS \(angstrom\)\n((?:.*\n){3})', s)
        if cm: cell = np.array([[float(v) for v in l.split()] for l in cm.group(1).strip().splitlines()])
        pm = re.search(r'ATOMIC_POSITIONS \(crystal\)\n((?:\s*\S+\s+\S+\s+\S+\s+\S+.*\n)+)', s)
        if pm:
            rows = [l.split() for l in pm.group(1).strip().splitlines()]; syms = [r[0] for r in rows]; frac = np.array([[float(v) for v in r[1:4]] for r in rows])
    wall = re.search(r'PWSCF\s+:\s+(.*?)\s+CPU\s+(.*?)\s+WALL', txt)
    converged = 'convergence has been achieved' in txt and ('JOB DONE' in txt)
    nk = re.search(r'number of k points=\s+(\d+)', txt); nb = re.search(r'number of Kohn-Sham states=\s+(\d+)', txt); npw = re.search(r'Parallelization info.*?\n.*?\n\s+sum\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', txt, flags=re.S)
    ram = re.search(r'Estimated total dynamical RAM >\s+([\d.]+)\s+(\w+)', txt)
    bfgs = 'bfgs converged' in txt
    return dict(energies_Ry=E, E_final_Ry=E[-1] if E else None, forces_eV_A=forces, F_final=forces[-1] if forces else None, stress_kbar=stress, stress_final=stress[-1] if stress else None,
                final_cell_A=cell, final_frac=frac, final_symbols=syms, converged=converged, bfgs_converged=bfgs, wall=wall.group(2) if wall else None,
                n_kpoints=int(nk.group(1)) if nk else None, n_bands=int(nb.group(1)) if nb else None, ram=(ram.group(1) + ' ' + ram.group(2)) if ram else None, job_done='JOB DONE' in txt)


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
