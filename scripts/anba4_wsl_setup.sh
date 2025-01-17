#!/bin/bash


# install conda in wsl 
install_conda() {
    CONDA_DIR=$1
    if [ ! -d "$CONDA_DIR" ]; then
        mkdir -p "$CONDA_DIR"
        wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O "$CONDA_DIR/miniforge.sh"
        bash "$CONDA_DIR/miniforge.sh" -b -u -p "$CONDA_DIR"
        rm "$CONDA_DIR/miniforge.sh"
    else
        echo "Miniforge is already installed."
    fi

    source "$CONDA_DIR/bin/activate"
    conda init --all
    conda install pip
}


# Check if conda is on PATH and works
check_conda() {
    if command -v conda &> /dev/null; then
        echo "Conda is installed and on PATH."
    else
        echo "Conda is not installed or not on PATH."
        exit 1
    fi
}


setup_anba4_env() {
    ANBADIR=$1
    # Check if the anba4-env environment exists
    if conda env list | grep -q 'anba4-env'; then
        echo "The anba4-env environment already exists."
    else
        echo "Creating the anba4-env environment."
        conda create -n anba4-env -y fenics=2019.1.0 python=3.9
    fi

    # Activate the anba4-env environment
    conda init
    conda activate anba4-env

    if [ ! -d "$ANBADIR/anba4" ]; then
        cd "$ANBADIR"
        git clone https://github.com/ANBA4/anba4.git
    else
        cd "$ANBADIR/anba4"
        git pull 
    fi

    cd "$ANBADIR/anba4"
    pip install -e .
    conda install pyyaml pyvista
}




setup_b3p_env() {
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
}

install_b3p() {
    WORKDIR=$1
    # Check if the b3p directory exists under workdir
    if [ ! -d "$WORKDIR/b3p" ]; then
        echo "Cloning the b3p repository."
        git clone -b doc https://github.com/wr1/b3p.git "$WORKDIR/b3p"
        cd "$WORKDIR/b3p"
    else
        echo "The b3p repository already exists. Pulling the latest changes."
        cd "$WORKDIR/b3p"
        git pull 
    fi
    pip install -e .
}


run_b3p_anba_example() {
    WORKDIR=$1
    # navigate to b3p dir /examples
    cd $WORKDIR/b3p/examples
    # run a quick check to see if the installation was successful
    b3p build blade_test.yml
    b3p 2d blade_test.yml
}


main() {
    WORKDIR=~/b3p_wsl

    # Check if the workdir directory exists, if not, create it
    if [ ! -d "$WORKDIR" ]; then
        mkdir -p "$WORKDIR"
    fi

    check_conda || install_conda "$HOME/miniforge3"

    setup_anba4_env $WORKDIR

    setup_b3p_env

    install_b3p $WORKDIR

    run_b3p_anba_example $WORKDIR
}

main