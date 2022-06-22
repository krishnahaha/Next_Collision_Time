#Script file to compile, pass parameters and run the simulation for a hard sphere event driven molecular dynamics simulation.

rm edmdMain.o
g++ -o3 edmdMain.cpp -o edmdMain.o

#for npa in 10000 25000 50000 100000
#for npa in 75000
#for npa in 100 200
#for npa in 10000 25000 50000 75000 100000
for npa in 100
#for npa in 100000
do
  #for ph in 0.0010 0.0100 0.1 0.1500 0.2000 0.2500 0.3000
  #for ph in 0.0500
  #for ph in 0.0010 0.0050 0.0100 0.0500 0.1000 0.1500 0.2000 0.2500 0.3000
for ph in 0.0010
do
    cd ./functions
    rm parameter.par

    echo "npart ${npa}">>parameter.par
    echo "dim 3">>parameter.par
    echo "nColl 20000">>parameter.par
    echo "nInit 1">>parameter.par
    echo "phi ${ph}">>parameter.par
    echo "vInit 1">>parameter.par
    echo "nCollSave 1000">>parameter.par
    echo "nCollIg 10000">>parameter.par

    cd ../
    #valgrind --track-origins=yes --leak-check=full ./cellDivEdmd.o './functions/parameter.par'
    ./edmdMain.o './functions/parameter.par'
    #>>outputFile_Phi_${ph}_npart_${npa}_withoutO3Run2 &
    #sleep 10
  done
done
