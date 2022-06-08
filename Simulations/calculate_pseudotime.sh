echo "Write the name of the simulation (e.g. velosim or dyngen)"
read method

if [[ $method == "velosim" ]]; then
    echo "Write the type of experiment (bifurcation, dropout or nodinamics)"
    read experiment
    if [[ $experiment == "dropout" ]]; then
        echo "Write the level of dropout you generated"
        read dropout
        arr=( ./datasets/$experiment/$method/$dropout/*.rds )
        resultsPath=./results/$experiment/$method/$dropout
    elif [[ $experiment == "bifurcation" ]]; then
        arr=( ./datasets/$experiment/$method/*.rds )
        resultsPath=./results/$experiment/$method
    elif [[ $experiment == "nodinamics" ]]; then
        arr=( ./datasets/no-dinamics/$method/*.rds )
        resultsPath=./results/no-dinamics/$method
    else 
        echo "You have chosen neither dropout, bifurcation or nodinamics"
        echo "Please, try again"
        exit
    fi
elif [[ $method == "dyngen" ]]; then
    echo "Write the type of experiment (linear or bifurcation)"
    read experiment
    if [[ $experiment == "linear" ]]; then
        arr=( ./datasets/$experiment/$method/*.rds )
        resultsPath=./results/$experiment/$method
    elif [[ $experiment == "bifurcation" ]]; then
        arr=( ./datasets/$experiment/$method/*.rds )
        resultsPath=./results/$experiment/$method
    else
        echo "You have chosen neither linear nor bifurcation"
        echo "Please, try again"
        exit
    fi
else 
    echo "You have chosen neither Velosim nor Dyngen"
    echo "Please, try again"
    exit
fi

echo "slingshot,paga" > $resultsPath/ti_pseudotime_cors.csv
echo "velocyto,scvelo" > $resultsPath/velos_pseudotime_cors.csv

n=1
for i in "${arr[@]}"
do
    Rscript ./code/preprocess_and_save_h5.r $i pseudotimeTrue
    python ./code/plots_velocity.py $i pseudotimeTrue

    n=$((n+1))
    if [[ $n -eq 50 ]]; then
        break
    fi
done

Rscript ./code/plot_pseudotime.r $resultsPath
