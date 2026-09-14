# rcm-spr-test4b-compute

Frozen Quantum ESPRESSO (7.5, PBEsol, PSlibrary PAW) input package and GitHub Actions workflows for a GaAs/AlAs [001] coherent-interface finite-displacement phonon calculation (Level-C Test 4B compute stage). Contains only compute inputs, a numpy-only output validator and workflows; scientific post-processing happens elsewhere. Every job verifies the frozen `pw.in` and pseudopotential SHA-256 hashes before running and uploads `pw.out`, `pw.err`, `time.txt`, `result.json` as an artifact named after the job id.

Workflows: `test4b-qe-probe` (one job), `test4b-qe-displacements` (matrix from `runs/level_c_gaas_alas_test4b/queue/gh_matrix.json`), `test4b-qe-collect` (aggregate artifacts).
