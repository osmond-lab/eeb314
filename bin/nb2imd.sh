#!/bin/bash

usage() {
    cat <<EOM
    ==================
    |||   nb2imd   |||
    ==================   
    Description:
        Jupyter notebook to interactive markdown. Takes a jupyter notebook and 
        makes a markdown page with embedded, executable interactive python cells
    Usage:
        -i <input> : The path to the jupyter notebook
        -o <output> : The path to the output directory
    
EOM
    exit 0
}

# Provide documentation if no command-line arguments are passed
if [ $# -lt 1 ] ; then
    usage
fi

# Get command-line arguments
while getopts i:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

# Convert the jupyter notebook to markdown
jupyter nbconvert --to markdown --output-dir="${output}" ${input} --ExecutePreprocessor.enabled=False

# Get the name of the file from input
split=$(echo ${input} | rev | cut -d "/" -f 1 | rev)
split=$(echo ${split} | cut -d "." -f 1)

# Make the code cells interactive in the markdown document
echo ${split}
python3 ./add_interactivity.py -f "${output}/${split}.md" -o "${output}/${split}.md"
