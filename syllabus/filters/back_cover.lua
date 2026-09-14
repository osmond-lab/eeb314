-- Appends a final page to the printed PDF pointing to the live course
-- website, vertically centered with a real hyperlink. Print-only content,
-- not sourced from any of the site's markdown pages.

function Pandoc(doc)
  if FORMAT ~= 'latex' and FORMAT ~= 'beamer' then
    return doc
  end

  local page = pandoc.RawBlock('latex', table.concat({
    '\\clearpage',
    '\\vspace*{\\fill}',
    '\\begin{center}',
    '\\Large See \\url{https://osmond-lab.github.io/eeb314} for more info',
    '\\end{center}',
    '\\vspace*{\\fill}',
  }, '\n'))

  doc.blocks:insert(page)
  return doc
end
