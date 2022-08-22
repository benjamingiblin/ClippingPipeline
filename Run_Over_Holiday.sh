#!/bin/bash

sleep 48h
echo "I am sleeping for 48h"
scancel -u bengib
Tree_Correlation_Function/Clean_Workers.sh

# Check which maps failed for SS3.11 if any.
# note this is not checking fid cosmol - have to check that manually when you're back
SS=3.11
./Check_Missing_Mass_Maps.sh $SS
echo "I just set 2nd sweep of SS3.11 going. Now I am sleeping for 30h"
sleep 30h

# One last sweep to make sure everything was produced okay
scancel	-u bengib
Tree_Correlation_Function/Clean_Workers.sh
./Check_Missing_Mass_Maps.sh $SS
echo "I just set 3rd sweep of SS3.11 going. Now I am sleeping for 30h"
sleep 30h

# Now move onto next smoothing scale:
# (note this will not run fid cosmol - have to do that manually when back)
SS=84.85
scancel	-u bengib
Tree_Correlation_Function/Clean_Workers.sh
./Check_Missing_Mass_Maps.sh $SS
echo "I just set 1st sweep of SS84.85 going. Now I am sleeping for 30h"
sleep 30h


scancel	-u bengib
Tree_Correlation_Function/Clean_Workers.sh
./Check_Missing_Mass_Maps.sh $SS
echo "I just set 2nd sweep of SS84.85 going. Now I am sleeping for 30h"
sleep 30h


scancel	-u bengib
Tree_Correlation_Function/Clean_Workers.sh
./Check_Missing_Mass_Maps.sh $SS
echo "I just set 3rd sweep of SS84.85 going. Now I am sleeping for 30h"
sleep 30h


scancel -u bengib
Tree_Correlation_Function/Clean_Workers.sh



