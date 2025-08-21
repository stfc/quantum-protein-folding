from pathlib import Path
import os
import sys
import warnings
from dotenv import load_dotenv


def load_project_root():
    """Load the project root directory from the .env file. Call this function at the beginning of all scripts."""
    load_dotenv()
    project_root = os.getenv("PROJECT_ROOT")
    if not Path(project_root).exists():
        warnings.warn(
            f"Project root directory does not exist: {project_root}. Imports might not work. Please check your .env file."
        )
    else:
        sys.path.append(str(project_root))
