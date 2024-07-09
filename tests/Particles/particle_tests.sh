cwd=$(pwd)

cd ~/tum-mglet-base/tests/Particles/ParticleAdvection1
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/ParticleAdvection2
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/ParticleAdvection3
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/ParticleAdvectionDiffusion
./run.sh newtest ../../../build/src/mglet

cd ~/tum-mglet-base/tests/Particles/ParticleDiffusion
./run.sh newtest ../../../build/src/mglet

cd $cwd
