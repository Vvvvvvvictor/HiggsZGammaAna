for mode in "test" "val";
do
    for channel in "zero_to_one_jet" "two_jet";
    do
        for score_type in "score_t" "score";
        do
            for type in 'background' 'signal';
            do
                python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/compare_4fold.py --type $type --score_type $score_type --channel $channel --mode $mode
            done
        done
    done
done