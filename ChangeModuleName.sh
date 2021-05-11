#!/bin/bash

currentName='sloth'
newName='sloth'

grep -rl "${currentName}" ./ | xargs sed -i "s/${currentName}/${newName}/g"
