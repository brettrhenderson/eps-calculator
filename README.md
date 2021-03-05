# Eps Calculator

Calculate the relative permittivity of bulk solids from Quantum Espresso CPMD output.

## Setup
To get this working on a compute canada cluster, first create a virtual environment.

### 1. Virtual Environment Setup
```
mkdir ~/.virtualenvs
module load python/3.8.5
virtualenv --no-download ~/.virtualenvs/pybec
source ~/.virtualenvs/pybec/bin/activate
```

### 2. Download the repo
```git clone git@github.com:Paci-Group/eps-calculator.git```

### 3. Install Dependencies
With your virtualenv activated, 

```
pip install --no-index --upgrade pip
pip install -r requirements.txt
```

### 4. Make it Executable
The script `calculate_eps` will load the necessary modules (on niagara at least)
to run the python module and calculate the permittivity. Place `calculate_eps`
somewhere in your PATH, or add the current directory to your path. Then, edit
the path in `calculate_eps` to point to wherever you cloned the git repo.
For example, I cloned it inside a directory in my $HOME folder called `pyscripts`, 
so the path is 
> $HOME/pyscripts/calculate_eps/calculate.py

### 5. Run It
In one of your job directories, you should now be able to run 
```
calculate_eps 7_zero_field.out 9_efield_ion.out --clamped-ion 8_clamped_ion.out --plot
```
to get a printout of your constant. You can use globs (*) to include multiple relaxed
ion outputs as well.

## Usage
* Note: When using `--plot`, you must be using X-forwarding in your SSH session.

```
usage: calculate.py [-h] [--clamped-ion CLAMPED_ION] [--plot] [--efield EFIELD] [--eps-bulk EPS_BULK] [--eps-inf-bulk EPS_INF_BULK] [--inclusion-element INCLUSION_ELEMENT]
                    zero_field relaxed_ion [relaxed_ion ...]

Calculate the dielectric permittivity of a bulk material from Quantum Espresso polarization output files.

positional arguments:
  zero_field            The file containing the zero field polarization output.
  relaxed_ion           The file or files containing the relaxed ion polarization output.

optional arguments:
  -h, --help            show this help message and exit
  --clamped-ion CLAMPED_ION, -c CLAMPED_ION
                        The file containing the clamped ion polarization output.
  --plot, -p            Plot the extracted polarization curve to check it.
  --efield EFIELD       Electric field in au, default=0.001.
  --eps-bulk EPS_BULK   Rel Permittivity of the bulk, for calculating enhancement, default=9.26.
  --eps-inf-bulk EPS_INF_BULK
                        High frequency Rel Permittivity of the bulk, default=3.04
  --inclusion-element INCLUSION_ELEMENT
                        Inclusion element symbol, default=Ag
```