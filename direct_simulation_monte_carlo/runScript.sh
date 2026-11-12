make


#npa=1000000
npa=1000188 #Get the data for this number of particles.
lambdaMult=1
ntrials=1
maxTime=1 #This times the avgFFT
it=0
while((it < ntrials));
do 
    #for ph in 0.2000 0.2500 0.3000 0.3500 0.4000 0.45000
    #for ph in 0.0010 0.0100 0.1000 0.1500 0.2000 0.2500 0.3000
    for ph in 0.20
    do
        rm parameter.par

        echo "npart ${npa}">>parameter.par
        echo "phi ${ph}">>parameter.par
        echo "temp 1.0">>parameter.par
        echo "lambdaMult ${lambdaMult}">>parameter.par
        echo "maxTime ${maxTime}">>parameter.par
        echo "dtFrac 0.2">>parameter.par
        echo "wallTemp 1.0">>parameter.par
        echo "wallVel 0.2">>parameter.par
        echo "trial ${it}">>parameter.par

        ./dsmc parameter.par &
        sleep 2
    done
    ((it++))
    sleep 0
done
