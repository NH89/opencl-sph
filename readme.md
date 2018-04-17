OpenCL-Athena
==========
How to run on Cluster

first ssh
  ```shell
  ssh USERNAME@athena.pawsey.org.au

  ```
  
Clone   
   ```shell
  git clone https://github.com/NH89/opencl-sph
  ```
  
Load Intel modules for compatibility
  ```shell
  module load intel-opencl-sdk 
  ```
  
Create build and install folder, and navigate into that folder
   ```shell
  mkdir mybuild
  mkdir install
  cd mybuild
  ```

Run Cmake with hint
   ```shell
  cmake -DCMAKE_INSTALL_PREFIX=../install .. -DOpenCL_LIBRARY=/pawsey/opencl-sdk/7.0.0/opencl/SDK/lib64/libOpenCL.so
  ```
(Alternatively load cuda module)
   ```shell
  cmake -DCMAKE_INSTALL_PREFIX=../install .. 
  ```
make targets
  ```shell
  make
  make install 
  ```
Reserve resources of GPU
```shell
srun -A gpuhack02 -N 1 -p gpuq --reservation=gpu-hackathon --gres=gpu:2 --pty /bin/bash
```
Run executable
  ```shell
  ./install/bin/opencl_sph
  ```

Sit back and drink a coffee


OpenCL-SPH
==========

An implementation of smoothed particle hydrodynamics using OpenCL.

This code is under development and may be subject to heavy revision at any time.

The nearest-neighbour search (NNS) algorithm is derived from the method used for [Fluids v3. by Rama Hoetzlein](https://github.com/rchoetzlein/fluids3).

Installation
------------

* Clone this repository

  ```shell
  git clone https://github.com/tauroid/opencl-sph
  ```

* Edit the first line of [makefile](makefile) to point to your MATLAB installation

* `make`

You will need OpenCL headers and libraries.

This has only been tested on linux, with gcc.

Usage
-----

The current interface is through MATLAB - see [go.m](go.m) or [go_fluid.m](go_fluid.m).

These implement an elastic solid based on [Becker et al. 2009](http://cg.informatik.uni-freiburg.de/publications/2009_NP_corotatedSPH.pdf) and a simple pressure gradient fluid respectively.

Using particle configurations other than a cube at the moment requires manual specification in the `position` field of psdata.

Profiling information
==========

On Athena

1. checkout Athena branch
1. ```sh cd ${HOME}```
1. ```sh module load intel intel-opencl-sdk broadwell forge```
1. ```sh make-profiler-libraries```
1. ```srun -A gpuhack02 -N 1 -p gpuq --reservation=gpu-hackathon --gres=gpu:2 --pty /bin/bash```
1. cd to your opencl-sph directory
1. ```sh mkdir -p mybuild install```
1. ```sh cd mybuild```
1. ```sh ../profiling/build_for_map.sh```
1. ```sh make```
1. ```sh make install```
1. ```sh cd ../profiling```
1. edit profile_code.sh to point to your mybuild/install/bin/opencl-sph location
1. ```sh ./profile_code.sh```
1. copy the .map file to your laptop
1. download map from https://developer.arm.com/products/software-development-tools/hpc/downloads/download-arm-forge to your laptop
1. open local map file and enjoy

