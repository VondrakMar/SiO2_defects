#!/bin/bash
PYTHON_SCRIPT="runMD.py" 
prepos="mvondrak_SiO2hookMay_"
get_data_python_path=/nexus/posix0/FHI-Theory/mvondrak/SiO2_defect_hookMay8/
for hook_folder in ${prepos}*/
do
    cd $hook_folder
    for run in run*/
    do
        cd $run
        i=$(grep -oP 'i=\K\d+' "$PYTHON_SCRIPT" | head -1)
        j=$(grep -oP 'j=\K\d+' "$PYTHON_SCRIPT" | head -1)
        k=$(grep -oP 'k=\K\d+' "$PYTHON_SCRIPT" | head -1)
        r0=$(grep -oP 'r0=\K[\d.]+' "$PYTHON_SCRIPT" | head -1)
        filename=$(grep -oP 'read\(\K"[^"]+"' "$PYTHON_SCRIPT" | tr -d '"')
        {
            echo "$i"
            echo "$j"
            echo "$k"
            echo "$r0"
            echo "$filename"
        } > data_here.data
        pwd
        cp ${get_data_python_path}get_data.py .
        python get_data.py
        cd ..
    done
    cd ..
done
