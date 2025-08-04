for mode in "test" "val";
do
    for type in 'signal' 'background' ;
    do
        for channel in "two_jet"; # "zero_to_one_jet"
        do
            python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/compare_4fold.py --type $type --channel $channel --mode $mode --variable H_mass
            for score_type in "score" "score_t";
            do
                python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/compare_4fold.py --type $type --score_type $score_type --channel $channel --mode $mode --variable bdt_score
            done
        done
    done
done