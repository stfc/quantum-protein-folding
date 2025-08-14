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
