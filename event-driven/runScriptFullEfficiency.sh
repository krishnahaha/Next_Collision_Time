rm edmd
make
ntrials=2

maxtime=("100000" "20000" "10000" "7000" "4200" "1700" "800" "500" "180" "80" "20" "13" "6" "3.5")
sampleObsInterval=("100001" "20001" "10001" "7001" "4201" "1701" "801" "501" "181" "81" "21" "13.5" "6.5" "3.6")

i=0
it=0

for ph in 0.4500;
do  
    i=0
    for npa in 108 500 ## 864 1372 2048 5324 10976 19652 48668 97556 275684 389344 780448 1314036;
    do
        while((it < ntrials));
        do
            rm parameter.par

            fileName=snapshot_end_n${npa}_pf${ph}_maxTime1000.0.dat
            echo "npart ${npa}">>parameter.par
            echo "phi ${ph}">>parameter.par
            echo "temp 1.0">>parameter.par
            echo "maxtime ${maxtime[$i]}">>parameter.par
            echo "sampleObsInterval ${sampleObsInterval[$i]}">>parameter.par
            echo "eps 1">>parameter.par
            echo "posInitType 1">>parameter.par #0: Random, 1: FCC, 2: File
            echo "shellSize 1.5">>parameter.par
            echo "fileName ${fileName}">>parameter.par
            echo "useThermostat 0">>parameter.par
            echo "trial ${it}">>parameter.par
            ((it++))

            ./edmd parameter.par & # >>overlapCheck_npa${npa}_phi${ph}.dat &
            sleep 10
            ((it++))
        done
        ((i++))
        #sleep 10
    done
done
