from pathlib import Path

ROOT_DIR = Path(__file__).parent.parent.parent

# data

DATA_DIR = ROOT_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"

# data/raw
PDB_DIR = RAW_DATA_DIR / "pdb"

# data/processed

# pyrosetta data
PYROSETTA_DIR = PROCESSED_DATA_DIR / "pyrosetta"
PYROSETTA_ENERGY_FILES_DIR = PYROSETTA_DIR / "pyrosetta_energy_files"
PYROSETTA_ENERGY_IMMUTABLE_FILES_DIR = PYROSETTA_DIR / "pyrosetta_energy_files_immutable"    # files used for simulations in the paper

# exact data
EXACT_ENERGY_DATA_DIR = PROCESSED_DATA_DIR / "exact"

# qaoa data
QAOA_DIR = PROCESSED_DATA_DIR / "qaoa"
MPS_QAOA_FILE = QAOA_DIR / "mps-qaoa.csv"
SV_QAOA_FILE = QAOA_DIR / "sv-qaoa.csv"

# simulated annealing data
SIMULATED_ANNEALING_DIR = PROCESSED_DATA_DIR / "annealing"
SIMULATED_ANNEALING_FILE = SIMULATED_ANNEALING_DIR / "sa.csv"
# qaoa paperv1 data
QAOA_PAPERV1_DIR = PROCESSED_DATA_DIR / "qaoa_paperv1"
MPS_QAOA_PAPERV1_FILE = QAOA_PAPERV1_DIR / "mps-qaoa.csv"
SV_QAOA_PAPERV1_FILE = QAOA_PAPERV1_DIR / "sv-qaoa.csv"

# simulated annealing data
SIMULATED_ANNEALING_DIR = PROCESSED_DATA_DIR / "annealing"
SIMULATED_ANNEALING_FILE = SIMULATED_ANNEALING_DIR / "sa.csv"

# simulated annealing paperv1 data
SIMULATED_ANNEALING_PAPERV1_DIR = PROCESSED_DATA_DIR / "annealing_paperv1"
SIMULATED_ANNEALING_PAPERV1_FILE = SIMULATED_ANNEALING_PAPERV1_DIR / "sa.csv"
