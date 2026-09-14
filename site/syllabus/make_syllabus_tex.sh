# Regenerate the cover image (site favicon) as a PDF pdflatex can embed.
inkscape "$(readlink -f ../images/favicon.svg)" --export-type=pdf --export-filename="$(readlink -f img/favicon.pdf)"

pandoc -s \
  --lua-filter=filters/cover.lua \
  --lua-filter=filters/strip_examples.lua \
  --lua-filter=filters/strip_textbook_links.lua \
  --lua-filter=filters/admonitions.lua \
  --lua-filter=filters/back_cover.lua \
  --include-in-header=filters/admonition_preamble.tex \
  -o syllabus.tex \
  ../index.md instructors.md when_and_where.md course_structure.md textbook.md final_project.md

# mkdocs (pymdownx.tasklist) renders "- [x] ..." items as a checkmark in a
# box. Pandoc instead renders checked task items as a box with an X
# ($\boxtimes$), the opposite look from a completed item -- swap it for a
# boxed checkmark so the PDF matches the site. Unchecked items already
# render as a plain empty box ($\square$), which matches fine as-is.
sed -i 's/\\item\[\$\\boxtimes\$\]/\\item[$\\textcolor{taskgreen}{\\boxed{\\checkmark}}$]/g' syllabus.tex
