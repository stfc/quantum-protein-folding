# quantum-protein-folding

A repository to accompany the paper [Quantum Algorithm for Protein Side-Chain Optimisation: Comparing Quantum to Classical Methods](https://arxiv.org/abs/2507.19383) by Anastasia Agathangelou, Dilhan Manawadu, and Ivano Tavernelli.

This repository provides:

* Sample Python scripts implementing classical (simulated annealing) and quantum (Quantum Approximate Optimisation Algorithm (QAOA)) algorithms to perform side-chain optimisation in proteins.
* Results and data supporting the findings of the paper.

## Directory structure

```
.
├── data
│   ├── processed
│   │   ├── annealing
│   │   │   ├── sa-paper.csv
│   │   │   └── sa.csv
│   │   ├── exact
│   │   │   ├── all-ground-state-energies.csv
│   │   │   └── res-*-rot-*.json
│   │   ├── pyrosetta
│   │   │   ├── pyrosetta_energy_files
│   │   │   │   ├── *rot_*res_one_body_terms.csv
│   │   │   │   └── *rot_*res_two_body_terms.csv
│   │   │   └── pyrosetta_energy_files_immutable
│   │   │       ├── *rot_*res_one_body_terms.csv
│   │   │       └── *rot_*res_two_body_terms.csv
│   │   └── qaoa
│   │       ├── mps-qaoa-paper.csv
│   │       ├── mps-qaoa.csv
│   │       ├── sv-qaoa-paper.csv
│   │       └── sv-qaoa.csv
│   ├── raw
│   │   ├── pdb
│   │   │   └── *residue.pdb
│   │   └── README.md
│   └── README.md
├── pyproject.toml
├── quantum_protein_folding
│   ├── quantum_protein_folding
│   │   ├── __init.py__
│   │   ├── bitstrings.py
│   │   ├── config.py
│   │   ├── ising.py
│   │   └── qaoa.py
│   ├── setup.py
│   └── src
│       ├── diagonal_estimator.py
│       ├── README.md
│       └── sampling_vqe.py
├── README.md
├── requirements.txt
└── scripts
    ├── calculate_exact_energies.py
    ├── create_energy_files.py
    ├── path_setup.py
    ├── run_qaoa.py
    └── run_simulated_annealing.py
```
### Installation

To set up the project locally, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/stfc/quantum-protein-folding.git
   cd quantum-protein-folding
   ```
2. Create a virtual environment using Python 3.11:
   ```
   python3.11 -m venv venv-pf
   ```
3. Activate the virtual environment:
   -On Linux/macOS
   ```
   source venv-pf/bin/activate
   ```
   -On Windows
   ```
   venv-pf\Scripts\activate
   ```
4. Install the required dependencies
   ```
   pip install -r requirements.txt
   ```

PyRosetta is required to generate one- and two-body rotamer interaction energies. It is included in `requirements.txt`, but if the `pip` installation fails, you can install it manually:

```bash
pip install pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
```
Refer to [PyRosetta documentation](https://www.pyrosetta.org/downloads) for more information.

Alternatively, `data/processed/pyrosetta` contains several precomputed rotomer interaction terms.

After the installation, create a `.env` file and specify your `PROJECT_ROOT`. This allows execution of scripts within the `scripts` directory.

## Running simulations

### Calculating one and two body rotomer interaction energies using pyrosetta

To generate rotamer energy files, run:
```
python scripts/create_energy_files.py -p <pdb_filename> -r <num_rotamers> -i <starting_rotomer_index>
```
Please refer to `scripts/create_energy_files.py` for more details including information on command line arguments.

### Calculating exact energies of bitstrings for benchmarking

The exact energies can be found by brute-force search over all allowed bitstrings. To do so, run
```
python scripts/calculate_bitstring_energies.py -res <num_res> -rot <num_rot>
```
The energy distributions are recorded in `json` format, and saved to `data/processed/exact` directory. Please note that this script uses as input the energy files in `data/processed/pyrosetta/pyrosetta_energy_files_immutable` directory. This can be changed by editing `quantum_protein_folding/quantum_protein_folding/config.py`.

### Running QAOA simulations

Execute:
```
python scripts/run_qaoa.py -res <num_res> -rot <num_rot> -p <num_ansatz_layers> -alpha <CVaR_aggregation> -s <num_shots> -m <simulator>
```
where the `simulator` can be either `MPS` or `SV` for `matrix-product-state` and `statevector` simulations, respectively.

### Running simulated annealing

Execute:
```
python scripts/run_simulated_annealing.py -res <num_res> -rot <num_rot>
```

## Data

Please refer to the `README.md` in `data` directory for more details.

## Authors

- [@dilhanm](https://github.com/DilhanM)
- [@anastasiaangelo](https://github.com/anastasiaangelo)

> [!NOTE]
> If you use this repository or its data in your work, we kindly request you to cite our paper with the following `bibtex` handle ([CITATION.bib](./CITATION.bib)):
> ```bibtex
>@ARTICLE{agathangelou2025quantumalgorithmproteinsidechain,
>      title={Quantum Algorithm for Protein Side-Chain Optimisation: Comparing Quantum to Classical Methods}, 
>      author={Anastasia Agathangelou and Dilhan Manawadu and Ivano Tavernelli},
>      year={2025},
>      eprint={2507.19383},
>      archivePrefix={arXiv},
>      primaryClass={quant-ph},
>      url={https://arxiv.org/abs/2507.19383}, 
>}
> ```

## Version History

* v0.1 - Initial Release