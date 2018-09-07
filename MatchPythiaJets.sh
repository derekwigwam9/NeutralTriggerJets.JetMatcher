#!/bin/bash
# 'MatchPythiaJets.sh'
#
# Use this to run 'MatchPythiaJets.C' in batch mode


# input files
parFile="\"input/pp200py.resTestPar.et920pi0.r02rm1chrg.d6m9y2018.root\""
detFile="\"input/pp200py.resTestDet.et920pi0.r02rm1chrg.d6m9y2018.root\""

# output file
outFile="\"pp200py.resTestResponse.r02a005rm1chrg.dr02q015185root\""


# run script
root -b -l <<EOF
  .L MatchPythiaJets.C++
  MatchPythiaJets($parFile, $detFile, $outFile, true)
  .q
EOF
