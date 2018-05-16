#!/bin/bash
### script altering fcs filenames to go through R script in right order
### example of use: bash Rename.sh 180222/plate1/
### example of use: bash Rename.sh /Users/mvlkova/Documents/RecA_Promoter/FC/180222/plate1/

### check if desired directory exists, if not, print error message
if [ ! -d $1 ]; then
  exit "
    Absolute or relative path to analyse data in must be provided as the first argument!
    Full path example: /Users/mvlkova/Documents/RecA_Promoter/FC/180222/plate1/
    Relative path example: 180222/plate1/"
fi

### filter files to be altered
aim="$(echo $1*_*_*_*.fcs | grep -v "*")"

### loop through filtered files
for file in $aim
do
  end=${file#*\/*\/*_*_*_}      ### save last number & file ending
  fin="${file%_*_*.fcs}_$end"   ### create new name of the file
  mv $file $fin                 ### rename the file
  echo "$file ---> $fin"        ### print info about the process
done
