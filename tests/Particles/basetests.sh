cwd=$(pwd)

cd ~/tum-mglet-base/tests/Particles/Advection1
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/Advection2
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/Advection3
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/Diffusion1
./run.sh newtest ../../../build/src/mglet

cd $cwd
