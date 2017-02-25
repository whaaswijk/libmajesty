# libmajesty
Library for Logic Synthesis and Optimization


## Installation

```bash
# As root, install the packages needed

# If the include path to Anaconda cannot be found, use the CPATH environment variable

# Boost
conda install -c conda-forge boost=1.61.0
echo 'export BOOST_ROOT=~/anaconda3/envs/<your-environment>' >> ~/.bashrc

# Cirkit
git clone https://github.com/msoeken/cirkit
cd cirkit/addons
git clone https://bitbucket.org/whaaswijk/cirkit-addon-experimental.git
cd ..
mkdir build
cd build
cmake ..
ccmake ..
# Turn on the formal and experimental addons!! Then press c to configure and q to quit
make external -j48
make cirkit -j48

# Set BOOST_ROOT in order to build Majesty!
echo 'export BOOST_ROOT=~/anaconda3'

# Finally, build Majesty
cd libmajesty
mkdir build
cmake ..
make

```

## Enabling Cirkit & ABC Features 

In order to enable the features offered by Cirkit and ABC (e.g. exact
MIG synthesis or resyn2), make sure that the Cirkit and ABC binaries
are on the path. Note that the ABC binary is automatically built by
Majesty.

```
# You can find the ABC binary under the Majesty build folder
cd libmajesty/build/abc/src/abc-project
pwd
# Add the output of pwd to the path
export PATH=$PATH:$(PWD_OUTPUT)
