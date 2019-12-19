Running the following commands should allow you to use the acados MEX interace on a linux machine.

```
git clone https://github.com/acados/acados.git
cd acados
git submodule update --recursive --init
mkdir build
cd build
cmake .. -DACADOS_WITH_QPOASES=ON -DACADOS_EXAMPLES=ON -DACADOS_UNIT_TESTS=ON
make -j4
make install -j4
# set paths to use acados MEX interface
cd ../examples/acados_matlab_octave/getting_started
source env.sh
# start matlab from this terminal
matlab &
```