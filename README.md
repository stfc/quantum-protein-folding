# quantum-protein-folding

## Scripts

### Calculating one and two body rotomer interaction energies using pyrosetta

1. Create a `.env` file and specify your `PROJECT_ROOT`. This allows execution of scripts within the `scripts` directory.
2. Install PyRosetta within your virtual environment.
3. Run the script to generate rotamer energy files:

   ```
   python scripts/create_energy_files.py -p <pdb_filename> -r <num_rotamers> -i <starting_rotomer_index>
   ```
   Please refer to `scripts/create_energy_files.py` for more details including information on command line arguments.