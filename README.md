# libmajesty
Library for Logic Synthesis and Optimization

TODO: remove _CRT_SECURE_NO_WARNINGS preprocessor flag and rewrite unsafe code


## Installation

```bash
# As root, install the packages needed


# Boost
conda install -c conda-forge boost=1.61.0
echo 'export BOOST_ROOT=~/anaconda3/envs/<your-environment>' >> ~/.bashrc

git clone https://github.com/redis/hiredis
cd hiredis
make
cd ..

hg clone https://bitbucket.org/alanmi/abc
cd abc
ABC_USE_PIC=1 make libabc.a -j24
cd ..

echo 'export HIREDIS_HOME=~/hiredis' >> ~/.bashrc
echo 'export LIBABC_HOME=~/abc' >> ~/.bashrc
```