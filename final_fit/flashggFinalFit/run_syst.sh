baseDir=`pwd`

cd $baseDir/Trees2WS/

# Set the input directory
inputDir="/eos/user/j/jiehan/root/output_cor_syst"

# Traverse all .root files and generate commands
find "$inputDir" -type f -name "*.root" | while read -r file; do
    # Extract the production mode and year
    proc_name=$(basename "$(dirname "$(dirname "$file")")")
    year_num=$(basename "$(dirname "$file")")

    # Generate and execute the command
    cmd="python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile $file --productionMode $proc_name --year $year_num"
    echo "$cmd"
    $cmd
done

# python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile /eos/user/j/jiehan/root/output_cor_syst --productionMode ggH --year 2017

# cd $baseDir/Signal/
# python3 RunSignalScripts.py --inputConfig config_ggH_2017.py --mode fTest
# python3 RunSignalScripts.py --inputConfig config_ggH_2017.py --mode calcPhotonSyst