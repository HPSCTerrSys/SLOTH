#!/bin/bash

currentName='catchyNAME'
newName='SOMENEWMDULENAME'

grep -rl "${currentName}" ./ | xargs sed -i "s/${currentName}/${newName}/g"
