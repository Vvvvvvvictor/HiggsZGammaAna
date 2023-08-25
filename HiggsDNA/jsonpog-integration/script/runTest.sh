#!/bin/bash
NOCHANGED_FILES=""
ADDED_FILES=""
CHANGED_FILES=""
BROKEN_FILES=""

pip install correctionlib==v2.0.0

#for i in $(find POG | grep "\.gz"); do # Run gunzip of .gz files
#    echo gunzip $i
#    gunzip $i
#done

for i in $(find POG | grep "\." | grep -v README.md); do # Whitespace-safe but not recursive.
    STATUS="$(correction validate --version 2 $i; echo $?)"
    if [[ ${STATUS: -1} -ne 0 ]]; then
        echo
        echo "######### ERROR in "$i" #########"
        correction validate --version 2 $i
        echo "#################################################################"
        echo
        BROKEN_FILES=$BROKEN_FILES"\n"$i
    else
        DIFF="$(cmp --silent cms-nanoAOD-repo/$i $i; echo $?)"
        if [[ $DIFF -ne 0 ]]; then
            if test -f "cms-nanoAOD-repo/$i"; then
                echo "There are changes in "$i" wrt cms-nanoAOD/jsonpog-integration.git. "
                echo "-------------- summary - original version -----------------------------------"
                correction summary cms-nanoAOD-repo/$i
                echo "-------------- summary - new version ----------------------------------------"
                correction summary $i
                echo "-------------- differences in summary ---------------------------------------"
                correction summary cms-nanoAOD-repo/$i | grep -v "Corrections in file" > tmp1.txt 
                correction summary $i | grep -v "Corrections in file" > tmp2.txt
                git diff --no-index tmp1.txt tmp2.txt
                echo "-------------- differences in file ----------------------------------------"
                git diff --no-index cms-nanoAOD-repo/$i $i
                echo "----------------------------------------------------------------------------"
                CHANGED_FILES=$CHANGED_FILES"\n"$i
            else
                echo "New file found in "$i
                echo "-------------- summary of new file -----------------------------------"
                correction summary $i
                echo "----------------------------------------------------------------------------"
                ADDED_FILES=$ADDED_FILES"\n"$i
            fi
        else
            echo "No changes in "$i" wrt cms-nanoAOD/jsonpog-integration.git. "
            NOCHANGED_FILES=$NOCHANGED_FILES"\n"$i
        fi
    fi
done

echo
echo -e "Files checked:"$NOCHANGED_FILES$CHANGED_FILES$ADDED_FILES$BROKEN_FILES
echo
echo -e "Files changed (tests passed):"$CHANGED_FILES
echo
echo -e "Files added (tests passed):"$ADDED_FILES
echo
if [[ ${#BROKEN_FILES} -ne 0 ]]; then
    echo -e "Broken files:"$BROKEN_FILES
    echo
    exit -1
else
    echo -e "No broken files."
fi
echo "Done."
