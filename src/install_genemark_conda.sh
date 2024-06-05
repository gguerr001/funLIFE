#!/bin/bash


# Define variables

GENEMARK_PACKAGE="src/gmes_linux_64_4.tar.gz"
INSTALL_DIR="$CONDA_PREFIX/genemark"

# Create installation directory
mkdir -p $INSTALL_DIR

# Download GeneMark-ES
wget $GENEMARK_URL -O $GENEMARK_PACKAGE

# Extract the package
tar -xzvf $GENEMARK_PACKAGE -C $INSTALL_DIR --strip-components=1

# Set environment variables
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo "export PATH=\$PATH:$INSTALL_DIR" > $CONDA_PREFIX/etc/conda/activate.d/genemark.sh
echo "export GENEMARK_PATH=$INSTALL_DIR" >> $CONDA_PREFIX/etc/conda/activate.d/genemark.sh

# Apply environment variables immediately for the current session
export PATH="$PATH:$INSTALL_DIR"
export GENEMARK_PATH="$INSTALL_DIR"


##Install Mutex
cpan App::cpanminus
cpanm MCE::Mutex


##Transfer key
#gunzip src/gm_key_64.gz
cp src/gm_key_64 ~/.gm_key

###AFTER RUNNING THIS
#Go to location where GeneMark is installed and run
original_dir=$(pwd)
cd "$CONDA_PREFIX/genemark" || exit
perl change_path_in_perl_scripts.pl "/usr/bin/env perl"


# Return to the original directory (if needed)
cd "$original_dir"




