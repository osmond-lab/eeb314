#!/bin/bash

# convert notebooks to markdown
for i in {1..9}; do ./nb2imd.sh -i ../notebooks/lectures/lecture-0$i.ipynb -o ../docs/lectures/; done
for i in {10..20}; do ./nb2imd.sh -i ../notebooks/lectures/lecture-$i.ipynb -o ../docs/lectures/; done
