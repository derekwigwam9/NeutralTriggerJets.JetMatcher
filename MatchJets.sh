#!/bin/bash
# 'MatchJets.sh'
#
# Use this to run 'MatchJets.C' in batch mode

# filepaths
gPath="\"input/geant/pp200r12pt9g.r04rm1full.d6m1y2018.root\""
uPath="\"input/mudst/pp200r12pt9u.r04rm1hc50full.d6m1y2018.root\""
oPath="\"pp200r12pt9.r04rm1hc50full.r04q15.d7m1y2018.root\""

root -b -l <<EOF
.L MatchJets.C++
MatchJets($gPath, $uPath, $oPath, true)
.q
EOF
