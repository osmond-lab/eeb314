# EEB314

This repository stores content for the Mathematical Modelling in Ecology and Evolution course (originally EEB430, now EEB314, at the University of Toronto) website. The website was created by @tomouellette -- many thanks!

## Building lectures

### Step 1

Write your lecture content into a jupyter notebook and save it into the `notebooks/lectures/` folder as `lecture-X.ipynb` (e.g. `lecture-01.ipynb`). If you have any images associated with this notebook, save them in a folder within `notebooks/lectures/` called `lecture-X-img/`.

### Step 2

To build your interactive book, now enter the `bin/` folder (i.e. `cd bin` from the root), and call the `nb2imd.sh` script as follows:

```bash
# Converts your notebook into an interactive markdown doc and moves it to correct folder
./nb2imd.sh -i ../notebooks/lectures/lecture-X.ipynb -o ../docs/lectures/

# If you have any images, move them to the correct folder
cp -r ../notebooks/lectures/lecture-X-img/ ../docs/lectures/lecture-X-img/
```

where `lecture-X` is the name of the lecture.

### Step 3

Open the `mkdocs.yml` file and go down to the line that says `nav`. The structure of the `nav` section determines the layout of the website. If this is a new lecture named `lecture-25.ipynb` then add it do the `nav` section as denoted below:

```yaml
nav:
  - Overview: index.md
  - Lectures:
    - Lecture 01: lectures/lecture-1.md
    - Lecture 25: lectures/lecture-25.md
```

Note that `.ipynb` is replaced with `.md` in the `nav` section.

## Building labs

### Step 1

The labs can be made using jupyter notebooks as before. 

### Step 2

When complete, they can be saved in a Google Drive folder.
 
### Step 3

To make them accessible to users, simply right click on the jupyter notebook in Google Drive and press `Open in Colaboratory`. Once in Google CoLab, click `Share` in the top right corner. Change general access to `Anyone with the link` and then copy the link in bottom right corner.

### Step 4

Take the copied/shareable link, and add it to the markdown table in `docs/labs/schedule.md`.

### Step 5

The labs can now be edited and saved directly in Google CoLab. Careful not to leave all the answers there!

## Building the site

### Step 1

You can view a local version of the website by calling `mkdocs serve` in the root of this folder (i.e. `mkdocs.yml` must be at the same level). Note you may need to install `mkdocs` and dependencies first with

```bash
pip install mkdocs
pip install mkdocs-material
```

### Step 2

If everything looks good locally, you can build the site calling `mkdocs build`. A new folder called `site` will be built that has everything you need for hosting the live website. You can now push to publish via GitHub Pages (see https://squidfunk.github.io/mkdocs-material/publishing-your-site/#with-github-actions for more info).

## More info

The website was built using [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/). Therefore, all the functionalitites described there can be added here if desired. Second, the executable cells were integrated using [Thebe](https://github.com/executablebooks/thebe). If any changes to the code blocks are desired, please look at their documentation for more information. Alternatively, you can directly modify the style of the code blocks by editing the `.css` stylesheets stored in `docs/stylesheets`.

### Site structure

The `nav` section controls the site structure. Adding or removing tags will add or remove the corresponding page. For example, if you wanted to only hide a section of the website you could go to the `nav` in `mkdocs.yml` and do:

```
nav:
  - Overview: index.md
  - Lectures:
    - Lecture 1: lectures/lecture-1.md
    - Lecture 2: lectures/lecture-2.md
    - Lecture 3: lectures/lecture-3.md
    # - Lecture 4: lectures/lecture-4.md
    # - Lecture 5: lectures/lecture-5.md
    # - Lecture 6: lectures/lecture-6.md
    # - Lecture 7: lectures/lecture-7.md
    # - Lecture 8: lectures/lecture-8.md
    # - Lecture 9: lectures/lecture-9.md
    # - Lecture 10: lectures/lecture-10.md
    # - Lecture 11: lectures/lecture-11.md
    # - Lecture 12: lectures/lecture-12.md
    # - Lecture 13: lectures/lecture-13.md
    # - Lecture 14: lectures/lecture-14.md
    # - Lecture 15: lectures/lecture-15.md
    # - Lecture 16: lectures/lecture-16.md
    # - Lecture 17: lectures/lecture-17.md
    # - Lecture 18: lectures/lecture-18.md
    # - Lecture 19: lectures/lecture-19.md
    # - Lecture 20: lectures/lecture-20.md
```

All the lines with `#` will be commented out and be removed from the site when you re-build. You can re-add them by uncommenting before re-building.
