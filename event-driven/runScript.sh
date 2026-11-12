rm edmd
make

npa=55296
ntrials=1
## maxtime=("10000" "1000" "100" "10" "10" "10" "10" "10" "10" "10")
## sampleObsInterval=("1" "0.1" "0.01" "0.001" "0.001" "0.001" "0.001" "0.001" "0.001" "0.001")

## maxtime=("1000" "1000" "1000" "1000" "1000" "1000" "1000" "1000" "1000" "1000")
## sampleObsInterval=("1" "1" "1" "1" "1" "1" "1" "1" "1" "1")


maxtime=("100" )
sampleObsInterval=("1")

it=0
i=0

#for ph in 0.2000 0.2500 0.3000 0.3500 0.4000 0.4500
## for ph in 0.0010 0.0100 0.1000 0.1500 0.2000 0.2500 0.3000 0.3500 0.4000 0.4500
## for ph in 0.2000 0.2500 0.3000 0.3500 0.4000 0.4500
for ph in 0.4500
do
    it=0
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
        echo "posInitType 2">>parameter.par #0: Random, 1: FCC, 2: File
        echo "shellSize 1.5">>parameter.par
        echo "fileName ${fileName}">>parameter.par
        echo "useThermostat 0">>parameter.par
        echo "trial ${it}">>parameter.par
        ((it++))

        ./edmd parameter.par >>screenPrint_npa${npa}_phi${ph}_maxTime${maxtime[$i]}.dat &
        ## ./edmd parameter.par >>viscosityAxx_npa${npa}_phi${ph}_maxTime${maxtime[$i]}.dat 
        sleep 10
    done
    ((i++))
done