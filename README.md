# libmajesty
Library for Logic Synthesis and Optimization

TODO: remove _CRT_SECURE_NO_WARNINGS preprocessor flag and rewrite unsafe code


## Installation

```bash
# As root, install the packages needed


# If the include path to Anaconda cannot be found, use the CPATH environment variable

# Boost
conda install -c conda-forge boost=1.61.0
echo 'export BOOST_ROOT=~/anaconda3/envs/<your-environment>' >> ~/.bashrc

git clone https://github.com/redis/hiredis
cd hiredis
make
cd ..

hg clone https://bitbucket.org/alanmi/abc
cd abc
ABC_USE_PIC=1 make libabc.a -j48
cd ..

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

echo 'export HIREDIS_HOME=~/hiredis' >> ~/.bashrc
echo 'export LIBABC_HOME=~/abc' >> ~/.bashrc

# Set BOOST_ROOT in order to build libmajesty!
echo 'export BOOST_ROOT=~/anaconda3'
```
