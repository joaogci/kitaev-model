#!/bin/bash
shopt -s extglob
 
if [ -e confout ]; then
  file_out=confout
  file_in=${file_out/out/in}
  echo "mv $file_out to $file_in"
  mv "$file_out" "$file_in"
fi
