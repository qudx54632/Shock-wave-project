# Quantum induced shock dynamics in gravitational collapse: insights from effective models and numerical frameworks

This repository contains the source code, documentation, and numerical results for a research project studying the formation and evolution of shock waves in spherically symmetric gravitational collapse, within an effective framework inspired by Loop Quantum Gravity (LQG).

The project explores how classical singularities are replaced by quantum-induced shell-crossing singularities and how these can be resolved via weak solutions, such as shock waves. A novel numerical scheme is developed to handle the nontrivial structures in the governing equations.

## Project Overview

- **Goal**: To model the evolution of shock waves in quantum-corrected black hole spacetimes and study the causal structure, horizon dynamics, and curvature properties.
- **Framework**: Effective Loop Quantum Gravity models using a generalized Painlevé–Gullstrand coordinate system.
- **Methods**: First-order partial differential equations, thin-shell junction conditions, and numerical shock-capturing schemes.

## Repository Structure

- `main.ipynb` — Main Jupyter notebook that runs the full simulation pipeline and generates figures.
- `.jl files` — Julia source scripts defining grid initialization, PDE solvers, flux updates, horizon detection, and characteristic tracking:
  - `initialize_grid.jl`, `get_PDE_solver.jl`, `flux_update.jl`, `update_state.jl`, etc.
  - `get_horizon.jl` and `characteristic_shock.jl` compute and track geometric features.
- `figures/` — Output plots and visualizations generated from the simulations.

## How to Run

To reproduce the results:

1. Open `main.ipynb` in Jupyter.
2. Make sure Julia is installed and the required packages are available.
3. The notebook will load and execute the Julia `.jl` files directly using embedded calls.
4. Figures will be saved into the `figures/` directory.

## Related Paper
If you are interested in the full theoretical background, please refer to the corresponding research paper: https://arxiv.org/pdf/...

## License

This project is made available under the [MIT License](LICENSE).

## Contact
For questions or collaborations, please contact:
	•	Name: Dongxue Qu
	•	Email: dqu@perimeterinstitute.ca
