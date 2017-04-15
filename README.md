# Majesty ♕
A Logic Synthesis and Optimization package


## Installation

```bash
# As root, install the packages needed
sudo apt-get install mercurial bison flex libreadline-gplv2-dev build-essential git g++ cmake cmake-curses-gui libgmp3-dev libxml2-dev libz-dev libncurses-dev

# Boost
wget --content-disposition https://downloads.sourceforge.net/project/boost/boost/1.63.0/boost_1_63_0.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fboost%2Ffiles%2Fboost%2F1.63.0%2F&ts=1492270618&use_mirror=netcologne
tar xf boost_1_63_0.tar.gz
cd boost_1_63_0
./bootstrap.sh
sudo ./b2 install

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

# Finally, build Majesty
git clone https://github.com/whaaswijk/libmajesty
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

```bash
# You can find the ABC binary under the Majesty build folder
cd libmajesty/build/abc/src/abc-project
pwd
# Add the output of pwd to the path
export PATH=$PATH:$(PWD_OUTPUT)
```