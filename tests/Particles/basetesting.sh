cwd=$(pwd)

echo ------------ ADVECTION 1 ----------
cd ~/tum-mglet-base/tests/Particles/Advection1
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION 2 ----------
cd ~/tum-mglet-base/tests/Particles/Advection2
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION 3 ----------
cd ~/tum-mglet-base/tests/Particles/Advection3
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION 4 ----------
cd ~/tum-mglet-base/tests/Particles/Advection4
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION 5 ----------
cd ~/tum-mglet-base/tests/Particles/Advection5
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION 6 ----------
cd ~/tum-mglet-base/tests/Particles/Advection6
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION 7 ----------
cd ~/tum-mglet-base/tests/Particles/Advection7
./run.sh newtest ../../../build/src/mglet

echo ------------ DIFFUSION 1 ----------
cd ~/tum-mglet-base/tests/Particles/Diffusion1
./run.sh newtest ../../../build/src/mglet

echo ------------ DIFFUSION 2 ----------
cd ~/tum-mglet-base/tests/Particles/Diffusion2
./run.sh newtest ../../../build/src/mglet

echo ------------ DIFFUSION 3 ----------
cd ~/tum-mglet-base/tests/Particles/Diffusion3
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION DIFFUSION 1 ----------
cd ~/tum-mglet-base/tests/Particles/AdvectionDiffusion1
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION DIFFUSION 2 ----------
cd ~/tum-mglet-base/tests/Particles/AdvectionDiffusion2
./run.sh newtest ../../../build/src/mglet

echo ------------ ADVECTION DIFFUSION 3 ----------
cd ~/tum-mglet-base/tests/Particles/AdvectionDiffusion3
./run.sh newtest ../../../build/src/mglet

cd $cwd
