#!/bin/bash

for D in *; do
    if [ -d "${D}" ]; then

    	if grep -Fxq "${D}" surfaces.txt
    	then
    		cd "${D}" # go to surface directory
    		blender -b -P ../anaglypy.py < ../options.txt
    		cd ..
    	# else
    		# do nothing
    	fi
    fi
done
