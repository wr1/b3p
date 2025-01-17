#!/bin/bash

# move to home
cd ~ 

# install conda in wsl 
if [ ! -d ~/miniforge3 ]; then
    mkdir -p ~/miniforge3
    wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O ~/miniforge3/miniforge.sh
    bash ~/miniforge3/miniforge.sh -b -u -p ~/miniforge3
    rm ~/miniforge3/miniforge.sh
else
    echo "Miniforge is already installed."
fi

source ~/miniforge3/bin/activate
conda init --all
conda install pip 


# Check if conda is on PATH and works
if command -v conda &> /dev/null; then
    echo "Conda is installed and on PATH."
else
    echo "Conda is not installed or not on PATH."
    exit 1
fi

# setup the default environment for anba and install anba in it 

# Check if the anba4-env environment exists
if conda env list | grep -q 'anba4-env'; then
    echo "The anba4-env environment already exists."
else
    echo "Creating the anba4-env environment."
    conda create -n anba4-env -y fenics=2019.1.0  python=3.9
fi

# Activate the anba4-env environment
conda init
conda activate anba4-env

cd ~
git clone https://github.com/ANBA4/anba4.git # (or git clone https://github.com/ANBA4/anba4.git)
cd anba4
pip install -e .
conda install pyyaml pyvista



# Check if the b3p environment exists with Python 3.11
if conda env list | grep -q 'b3p'; then
    echo "The b3p environment already exists. Activating it."
    conda activate b3p
else
    echo "Creating the b3p environment with Python 3.11."
    conda create -n b3p -y python=3.11
    conda activate b3p
    conda install pip
fi

# wget https://bootstrap.pypa.io/get-pip.py
conda install pip 


# Check if the projects directory exists, if not, create it
if [ ! -d ~/projects ]; then
    mkdir -p ~/projects
fi

# Check if the b3p directory exists under projects
if [ ! -d ~/projects/b3p ]; then
    echo "Cloning the b3p repository."
    git clone -b doc https://github.com/wr1/b3p.git ~/projects/b3p
    cd ~/projects/b3p
    git pull origin doc
else
    echo "The b3p repository already exists. Pulling the latest changes."
    cd ~/projects/b3p
    git pull origin doc
fi

pip install -e . 
pip install sympy 

# run a quick check to see if the installation was successful
cd examples
b3p build blade_test.yml
b3p 2d blade_test.yml




