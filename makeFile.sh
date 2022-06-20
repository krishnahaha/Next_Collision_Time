#Script file to compile, pass parameters and run the simulation for a hard sphere event driven molecular dynamics simulation.
#Edit the CMakeLists.txt file to set up the file dependencies needed. This change is needed whenever a compilation flags are modified.

cmake .
make
mv edmdMain Makefile cmake_install.cmake CMakeCache.txt ./build

#Creation of the parameter file.
for npa in 10
do
  for ph in 0.0010
  do
    cd ./functions
    rm parameter.par

    echo "npart ${npa}">>parameter.par
    echo "dim 3">>parameter.par
    echo "nColl 1000">>parameter.par
    echo "nInit 1">>parameter.par
    echo "phi ${ph}">>parameter.par
    echo "vInit 1">>parameter.par
    echo "nCollSave 1000">>parameter.par
    echo "nCollIg 10000">>parameter.par

    cd ../build
    #valgrind --track-origins=yes --leak-check=full ./edmdMain.o './functions/parameter.par'
    ./edmdMain '../functions/parameter.par'
    #>>outputFile_Phi_${ph}_npart_${npa}_O3Run1 &
    #sleep 10
    cd ../ #Come to the main directory.
  done
done
