cwd=$(pwd)

cd ~/tum-mglet-base/tests/Particles/Advection4
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/Advection5
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/Advection6
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/Diffusion2
./run.sh newtest ../../../build/src/mglet


cd $cwd
