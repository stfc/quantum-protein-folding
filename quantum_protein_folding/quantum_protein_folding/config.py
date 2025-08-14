from pathlib import Path

ROOT_DIR = Path(__file__).parent.parent.parent

# data

DATA_DIR = ROOT_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"

# data/raw
PDB_DIR = RAW_DATA_DIR / "pdb"

# pyrosetta data
PYROSETTA_ENERGY_FILES_DIR = RAW_DATA_DIR / "pyrosetta_energy_files"
PYROSETTA_ENERGY_FIXED_FILES_DIR = RAW_DATA_DIR / "pyrosetta_energy_files_immutable"    # files used for simulations in the paper