# for cat in {0..3}
# do
# # for sig in 0 1 2 5 10 20
# # do
# sig=0
# for chan in "zero_jet" "one_jet" "two_jet" "VH_ttH"
# do
# if [ -f "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.txt" ]; then
# rm "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.txt";
# touch "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.txt";
# fi
# # if [ -f "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.root" ]; then
# # rm "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.root";
# # fi
# for fun in "gauxpow1" "gauxpow3" "gauxpow5" "gauxexp1" "gauxexp3" "gauxexp5" "gauxlau1" "gauxlau2" "gauxlau3" "gauxlau4" "gauxlau5" "bern2" "bern3" "bern4" "bern5"
# do
# root -l -q 'SSTest/SSTest.cpp('$cat', '$sig', "'$chan'", "'$fun'")';
# done
# done
# # done
# done

for cat in {1..4}
do
# for sig in 0 1 2 5 10 20
# do
sig=0
# for chan in "zero_jet" "one_jet" "two_jet" "VH_ttH"
# do
chan="untagged"
if [ -f "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.txt" ]; then
rm "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.txt";
touch "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.txt";
fi
# if [ -f "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.root" ]; then
# rm "SSTest/outputs/"$chan"_"$cat"_"$sig"xsig.root";
# fi
for fun in "gauxpow1" "gauxpow3" "gauxpow5" "gauxexp1" "gauxexp3" "gauxexp5" "gauxlau1" "gauxlau2" "gauxlau3" "gauxlau4" "gauxlau5" "bern2" "bern3" "bern4" "bern5"
do
root -l -q 'SSTest/SSTest.cpp('$cat', '$sig', "'$chan'", "'$fun'")';
# done
done
# done
done