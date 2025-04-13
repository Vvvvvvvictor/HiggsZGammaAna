#!/bin/bash
TAKE_CARE_HOUR=9
# COMMAND="sh scripts/run_analysis_bkgmc_run3.sh"  
COMMAND="sh scripts/run_analysis_data_run3.sh"  
# COMMAND="sh scripts/run_analysis_signal_run3.sh"
SLEEP_TIME=5
RUN_TIME_MINUTES=15
RUN_TIME=$((RUN_TIME_MINUTES * 60))  # Convert minutes to seconds
MAX_TRIALS=$((TAKE_CARE_HOUR * 4))
TRIAL_COUNT=0

while [ $TRIAL_COUNT -lt $MAX_TRIALS ]; do
    clear
    echo "Trial $((TRIAL_COUNT + 1)) of $MAX_TRIALS starting at $(date)"
    $COMMAND
    PID=$!
    sleep $RUN_TIME
    echo "Attempting to suspend command (PID: $PID) at $(date)"
    kill -TSTP $PID 2>/dev/null
    
    # Check if process is still running (not terminated)
    if kill -0 $PID 2>/dev/null; then
        echo "Process $PID has been suspended at $(date)"
        echo "You can resume it later with 'fg' or 'bg' command"
    else
        echo "Process $PID no longer exists"
    fi
    
    TRIAL_COUNT=$((TRIAL_COUNT + 1))
    
    if [ $TRIAL_COUNT -eq $MAX_TRIALS ]; then
        echo "Reached maximum trials ($MAX_TRIALS). Exiting."
        break
    fi
    
    echo "Waiting $SLEEP_TIME seconds before next trial..."
    sleep $SLEEP_TIME
done

echo "Script completed at $(date)"