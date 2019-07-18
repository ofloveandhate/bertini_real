#!/bin/bash
# find . -maxdepth 1 -type d \( ! -name . \) -exec bash -c "cd '{}' && pwd" \;

for D in *; do
    if [ -d "${D}" ]; then

    	if grep -Fxq "${D}" surfaces.txt
    	then
    		cd "${D}" # go to surface directory
    		blender -b -P ~/bertini_real/python/bertini_real/anaglypy/animate.py < ../options.txt
    		cd ..
    	else
    		# do nothing
    	fi
    fi
done

