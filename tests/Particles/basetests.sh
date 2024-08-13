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

echo ------------ DIFFUSION 1 ----------
cd ~/tum-mglet-base/tests/Particles/Diffusion1
./run.sh newtest ../../../build/src/mglet

cd $cwd
