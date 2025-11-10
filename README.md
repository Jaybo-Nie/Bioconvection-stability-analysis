# Bioconvection Stability — Code & Data

This repository contains code and minimal data to reproduce the analyses for a study on the onset of bioconvection in suspensions of motile microorganisms.

## Repository layout
```
.
├── data/
│   ├── H11.csv                  # Results from bisection search for H = 2
│   └── Critical.csv             # Critical velocity at different layer depths
├── env/
│   ├── environment.yml          # Python/Conda environment
│   ├── Project.toml             # Julia environment file (pinned)
│   └── Manifest.toml            # Julia environment lockfile
├─ figures/                      # Generated figures (PDF)
├─ notebooks/
│   └─ stability_plot.ipynb      # Visualization of analytical and numerical thresholds
├─ simulations/                  # Julia FEM simulations and utilities
│   └─ Simulations.jl            # Main numerical simulation script
├─ LICENSE                       # MIT license
└─ README.md                     # this file
```


## Quick start (Julia)
```bash
# From the repository root
julia --project=env -e 'using Pkg; Pkg.instantiate()'

# Run the main simulation script
julia --project=env simulations/Simulations.jl

```
> Tip: The Julia environment is located under `env/`. If you move it to the project root, replace `--project=env` with `--project=`. accordingly.

## Reproducing figures
The workflow involves two main steps:

1. **Compute thresholds:** perform parameter sweeps and bisection searches; results are saved to the `data/` folder as CSV files.
2. **Generate plots** → process and visualize analytical and numerical thresholds using the Python notebook (`notebooks/stability_plot.ipynb`). Outputs are saved in `figures/`.

## Data
- `data/H11.csv` stores all results from the bisection search for depth H=2.
- `data/Critical.csv` contains the numerical critical velocities across multiple layer depths.

## Environment
- **Julia**: use the pinned environment defined in `env/Project.toml` and `env/Manifest.toml`.
- **Python**: use the Conda environment specified in `env/environment.yml`.
```bash
conda env create -f env/environment.yml

conda activate bioconvection

```
