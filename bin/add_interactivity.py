#!/bin/env python3

"""
Takes a jupyter notebook that has been converted to markdown and replaces
all code blocks with embedded, executable interactive cells using Thebe.
"""

import string
import argparse
import numpy as np

def arg_parse():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--file', type=str, help='The path to the file to process')
	parser.add_argument('-o', '--output', type=str, help='Where should the file be saved; the output directory and name of file')	
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	# Collect arguments
	file, output = vars(arg_parse()).values()

	# Load and read markdown file
	with open(file) as f:
		lines = f.readlines()

	# Setup header script to enable interactivity
	header_script = [
		'<script type="text/x-thebe-config">\n',
		'  {\n',
		'      requestKernel: true,\n',
		'      mountActivateWidget: true,\n',
		'      mountStatusWidget: true,\n',
		'      binderOptions: {\n',
		'      repo: "mmosmond/executable-cells",\n',
		'      ref: "main",\n',
		'      binderUrl: "https://gke.mybinder.org",\n'
		'      },\n',
		'  }\n',
		'</script>\n',
		'<script src="https://unpkg.com/thebe@latest/lib/index.js"></script>\n\n'
	]

	# Add header to lines
	lines = header_script + lines

	# Setup kernel for running code
	kernel = [
		'\n<hr style="margin-bottom: 0em;">\n',
		'<center>\n',
		'<div class="inrow">\n',
		'	Run notes interactively?\n',
		'	<div style="float: left;" class="thebe-activate"></div>\n',
		'	<div class="thebe-status"></div>\n',
		'</div>\n',
		'</center>\n',
		'<hr style="margin-top: 0em;">\n'
	]

	# Place kernel after header
	for i, n in enumerate(lines):
		if n[0] == '#':
			break

	# Add kernel after header
	lines = lines[0:(i+1)] + kernel + lines[(i+1):]

	# Find the code blocks in markdown file
	lines = np.array(lines)
	code_block_upper = np.where(lines == '```python\n')[0]
	code_block_lower = np.where(lines == '```\n')[0]

	# Make all code blocks interactive
	lines[code_block_upper] = '<pre data-executable="true" data-language="python">\n'
	lines[code_block_lower] = '</pre>\n'

	# Save updated markdown file
	with open(output, 'w') as f:
	    for line in list(lines):
	        f.write(line)



