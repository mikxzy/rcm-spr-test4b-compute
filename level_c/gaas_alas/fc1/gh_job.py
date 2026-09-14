"""Run ONE frozen Test 4B Quantum ESPRESSO job inside a GitHub Actions runner (or any Linux box with pw.x on PATH).

  python level_c/gaas_alas/fc1/gh_job.py --job bulkxv_GaAs_fd_000 --np 4 [--out artifacts/bulkxv_GaAs_fd_000]

Immutable job package = queue/jobs.json entry (job id, directory, calc type, stage, atom count, frozen pw.in sha256)
+ the frozen pw.in on disk + the pseudopotential files (sha256 recorded in TEST4B_PREREGISTRATION.json).
This script: verifies every hash BEFORE running, records the runner/software environment, runs pw.x under
/usr/bin/time -v with a memory sampler, then validates the output (JOB DONE, SCF/BFGS convergence, full force set, finite
values, stress) and writes result.json. It never edits pw.in. Requires only numpy (parser) — no phonopy.
Exit code: 0 = DONE (valid), 2 = QE ran but output invalid, 3 = QE non-zero exit, 4 = precondition (hash) failure."""
from __future__ import annotations
import argparse, hashlib, json, os, platform, shutil, subprocess, sys, time, threading
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT))
from level_c.gaas_alas.fc1.qe import parse_pwout                       # numpy-only parser (identical to the local one)

OUT4B = ROOT / 'runs/level_c_gaas_alas_test4b'; DFT = OUT4B / 'dft'; Q = OUT4B / 'queue'


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def env_record():
    def run(cmd):
        try: return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120).stdout.strip()
        except Exception as e: return f'ERR {e}'
    return dict(hostname=platform.node(), platform=platform.platform(), nproc=run('nproc'), cpu=run("lscpu | grep -E 'Model name|^CPU\\(s\\)|Thread|Core' | tr -s ' '"), meminfo=run("grep -E 'MemTotal|MemAvailable' /proc/meminfo"),
                pw_x=run('which pw.x'), pw_version=run("pw.x -h 2>&1 | head -2 | tail -1 || true"), mpirun=run('mpirun --version | head -1'), conda_packages=run("mamba list 2>/dev/null | grep -iE '^ *(qe|openmpi|libopenblas|libblas|liblapack|fftw|scalapack|elpa|python|numpy) '"),
                github=dict(runner=os.environ.get('RUNNER_NAME'), os=os.environ.get('RUNNER_OS'), arch=os.environ.get('RUNNER_ARCH'), image=os.environ.get('ImageOS'), run_id=os.environ.get('GITHUB_RUN_ID'), job=os.environ.get('GITHUB_JOB'), sha=os.environ.get('GITHUB_SHA'), workflow=os.environ.get('GITHUB_WORKFLOW')))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--job', required=True); ap.add_argument('--np', type=int, default=4); ap.add_argument('--out'); ap.add_argument('--max-seconds', type=int, default=0)
    ap.add_argument('--nk', type=int, default=1, help='QE k-point pools (-nk): pure MPI work distribution, no effect on the converged result')
    a = ap.parse_args(); jobs = json.load(open(Q / 'jobs.json')); job = jobs[a.job]; d = DFT / job['dir']; out = Path(a.out) if a.out else d / 'artifact'; out.mkdir(parents=True, exist_ok=True)
    pre = json.load(open(OUT4B / 'TEST4B_PREREGISTRATION.json')); res = dict(job=job, started=time.strftime('%Y-%m-%dT%H:%M:%S'), np=a.np, nk=a.nk, environment=env_record())
    # ---- preconditions: frozen input + pseudopotential hashes
    h = sha(d / 'pw.in'); res['pwin_sha256_observed'] = h; res['pseudopotentials'] = {}
    ok = h == job['pwin_sha256']
    for el, spec in pre['dft']['pseudopotentials'].items():
        hp = sha(DFT / 'pseudo' / spec['file']); res['pseudopotentials'][el] = dict(file=spec['file'], expected=spec['sha256'], observed=hp, ok=hp == spec['sha256']); ok = ok and hp == spec['sha256']
    if job['stage'] == 'fd':
        frozen = json.load(open(OUT4B / 'FC1_displacements_frozen.json')); i = int(job['id'].split('_')[-1]); rec = frozen['displacements'][i]
        res['displacement'] = rec; res['relaxed_structure_sha256'] = frozen['relaxed_structure_sha256']; ok = ok and rec['pwin_sha256'] == job['pwin_sha256']
    if not Path('/work/pseudo').exists(): res['note'] = '/work/pseudo missing (pw.in uses pseudo_dir=/work/pseudo)'; ok = False
    res['preconditions_ok'] = ok
    if not ok:
        res['status'] = 'FAILED'; res['reason'] = 'precondition failure (hash/paths)'; json.dump(res, open(out / 'result.json', 'w'), indent=1, default=str); print(json.dumps(res, indent=1, default=str)); sys.exit(4)
    # ---- run pw.x with a memory sampler
    for f in ['pw.out', 'pw.err', 'time.txt']:
        if (d / f).exists(): (d / f).unlink()
    samples = []; stop = threading.Event()
    def sampler():
        while not stop.is_set():
            try:
                mi = dict(l.split(':')[0:1] + [l.split()[1]] for l in open('/proc/meminfo') if l.startswith(('MemTotal', 'MemAvailable')))
                samples.append((round(time.time() - t0), int(mi['MemTotal']) - int(mi['MemAvailable'])))
            except Exception: pass
            stop.wait(10)
    gtime = next((t for t in [shutil.which('time'), '/usr/bin/time', '/opt/conda/bin/time'] if t and Path(t).exists()), None)
    timer = f'{gtime} -v -o time.txt ' if gtime else ''; res['gnu_time'] = gtime
    cmd = f'cd "{d}" && {timer}mpirun --oversubscribe --bind-to none -np {a.np} pw.x -nk {a.nk} -in pw.in > pw.out 2> pw.err'
    t0 = time.time(); th = threading.Thread(target=sampler, daemon=True); th.start(); p = subprocess.run(cmd, shell=True); stop.set(); th.join(timeout=1); wall = time.time() - t0
    subprocess.run('rm -rf /tmp/qe_*', shell=True)
    peak_rss_kb = None
    if (d / 'time.txt').exists():
        for l in open(d / 'time.txt'):
            if 'Maximum resident set size' in l: peak_rss_kb = int(l.split()[-1])
    res.update(exit_code=p.returncode, wall_s=round(wall, 1), peak_rss_per_rank_MB=round(peak_rss_kb / 1024, 1) if peak_rss_kb else None, peak_system_used_MB=round(max(s[1] for s in samples) / 1024, 1) if samples else None, mem_samples=samples[::6])
    # ---- validate
    r = parse_pwout(d / 'pw.out') if (d / 'pw.out').exists() else None; why = []
    if r is None: why.append('no pw.out')
    else:
        if not r['job_done']: why.append('no JOB DONE')
        if not r['converged']: why.append('SCF not converged')
        if r['F_final'] is None or (job['n_atoms'] and r['F_final'].shape[0] != job['n_atoms']): why.append('forces missing/incomplete')
        elif not bool(__import__('numpy').isfinite(r['F_final']).all()): why.append('non-finite forces')
        if r['stress_final'] is None: why.append('stress missing')
        if job['calc'] == 'vc-relax' and not (r['bfgs_converged'] and r['final_cell_A'] is not None): why.append('BFGS not converged / no final coordinates')
        if job['calc'] == 'vc-relax-chunk':
            txt = (d / 'pw.out').read_text(errors='ignore'); finished = r['bfgs_converged'] and r['final_cell_A'] is not None
            stopped = 'Maximum CPU time exceeded' in txt and 'ATOMIC_POSITIONS' in txt
            res['chunk'] = dict(finished=bool(finished), stopped_at_max_seconds=bool(stopped), n_bfgs_steps=len(r['energies_Ry']))
            if not (finished or stopped): why.append('chunk neither converged nor cleanly stopped at max_seconds')
        res.update(E_final_Ry=r['E_final_Ry'], n_scf_energies=len(r['energies_Ry']), n_kpoints=r['n_kpoints'], n_bands=r['n_bands'], qe_ram_estimate=r['ram'], qe_wall=r['wall'], max_force_eV_A=float(abs(r['F_final']).max()) if r['F_final'] is not None else None)
    res['integrity'] = dict(ok=not why, issues=why); res['status'] = 'DONE' if not why and p.returncode == 0 else 'FAILED'; res['reason'] = 'ok' if res['status'] == 'DONE' else '; '.join(why) or f'exit {p.returncode}'
    res['output_hashes'] = {f: sha(d / f) for f in ['pw.in', 'pw.out', 'pw.err', 'time.txt'] if (d / f).exists()}; res['ended'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    for f in ['pw.in', 'pw.out', 'pw.err', 'time.txt']:
        if (d / f).exists(): shutil.copy(d / f, out / f)
    json.dump(res, open(out / 'result.json', 'w'), indent=1, default=str); print(json.dumps({k: v for k, v in res.items() if k not in ('environment', 'mem_samples', 'job')}, indent=1, default=str))
    sys.exit(0 if res['status'] == 'DONE' else (3 if p.returncode != 0 else 2))


if __name__ == '__main__':
    main()
