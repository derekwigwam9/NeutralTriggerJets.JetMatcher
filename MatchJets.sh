#!/bin/bash
# 'MatchJets.sh'
#
# Use this to run 'MatchJets.C' in batch mode

# i/o filepaths
declare -a parFiles
declare -a detFiles
declare -a outFiles

# particle level input
parFiles[0]="\"../JetMaker/mc/pp200r9pt4rff.particle.r03rm1chrg.root\""
parFiles[1]="\"../JetMaker/mc/pp200r9pt5rff.particle.r03rm1chrg.root\""
parFiles[2]="\"../JetMaker/mc/pp200r9pt7rff.particle.r03rm1chrg.root\""
parFiles[3]="\"../JetMaker/mc/pp200r9pt9rff.particle.r03rm1chrg.root\""
parFiles[4]="\"../JetMaker/mc/pp200r9pt11rff.particle.r03rm1chrg.root\""
parFiles[5]="\"../JetMaker/mc/pp200r9pt15rff.particle.r03rm1chrg.root\""
parFiles[6]="\"../JetMaker/mc/pp200r9pt25rff.particle.r03rm1chrg.root\""
parFiles[7]="\"../JetMaker/mc/pp200r9pt35rff.particle.r03rm1chrg.root\""

# detector level input
detFiles[0]="\"../JetMaker/mudst/pp200r9pt4rff.et920vz55had.r03rm1chrg.root\""
detFiles[1]="\"../JetMaker/mudst/pp200r9pt5rff.et920vz55had.r03rm1chrg.root\""
detFiles[2]="\"../JetMaker/mudst/pp200r9pt7rff.et920vz55had.r03rm1chrg.root\""
detFiles[3]="\"../JetMaker/mudst/pp200r9pt9rff.et920vz55had.r03rm1chrg.root\""
detFiles[4]="\"../JetMaker/mudst/pp200r9pt11rff.et920vz55had.r03rm1chrg.root\""
detFiles[5]="\"../JetMaker/mudst/pp200r9pt15rff.et920vz55had.r03rm1chrg.root\""
detFiles[6]="\"../JetMaker/mudst/pp200r9pt25rff.et920vz55had.r03rm1chrg.root\""
detFiles[7]="\"../JetMaker/mudst/pp200r9pt35rff.et920vz55had.r03rm1chrg.root\""

# output
outFiles[0]="\"output/pp200r9pt4rff.matched.et920vz55had.r03rm1chrg.dr03q15.root\""
outFiles[1]="\"output/pp200r9pt5rff.matched.et920vz55had.r03rm1chrg.dr03q15.root\""
outFiles[2]="\"output/pp200r9pt7rff.matched.et920vz55had.r03rm1chrg.dr03q15.root\""
outFiles[3]="\"output/pp200r9pt9rff.matched.et920vz55had.r03rm1chrg.dr03q15.root\""
outFiles[4]="\"output/pp200r9pt11rff.matched.et920vz55had.r03rm1chrg.dr03q15.root\""
outFiles[5]="\"output/pp200r9pt15rff.matched.et920vz55had.r03rm1chrg.dr03q15.root\""
outFiles[6]="\"output/pp200r9pt25rff.matched.et920vz55had.r03rm1chrg.dr03q15.root\""
outFiles[7]="\"output/pp200r9pt35rff.matched.et920vz55had.r03rm1chrg.dr03q15.root\""


# loop over files
(( iFile=0 ))
for pFile in ${parFiles[@]}; do

# run script
root -b -l <<EOF
  .L MatchJets.C++
  MatchJets($pFile, ${detFiles[$iFile]}, ${outFiles[$iFile]}, true)
  .q
EOF
(( iFile++ ))

done


# delete arrays
unset parFiles
unset detFiles
unset outFiles
