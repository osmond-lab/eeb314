-- Builds the syllabus cover (page 1):
--   * a big centered cover image (the site favicon, pre-rendered to PDF in
--     img/favicon.pdf since pdflatex can't place SVGs directly)
--   * the course title -- originally the <h1> in index.md, which survives
--     pandoc's HTML handling only as a plain-text paragraph -- reset as a
--     big bold centered title
-- and forces a page break before the first top-level section so the
-- syllabus proper (Instructors, ...) starts on page 2. LaTeX/PDF output only.

local COVER_IMAGE = 'img/favicon.pdf'

function Pandoc(doc)
  if FORMAT ~= 'latex' and FORMAT ~= 'beamer' then
    return doc
  end

  local out = pandoc.Blocks({})
  local styled_title = false
  local inserted_break = false

  for _, b in ipairs(doc.blocks) do
    -- The <h1> in index.md survives pandoc's HTML handling as a bare
    -- Plain block (its surrounding <h1>/<table> tags are raw HTML that
    -- gets dropped by the LaTeX writer); an empty Plain also precedes it
    -- for the <img> tag (dropped the same way), so skip blank ones.
    if not styled_title and b.t == 'Plain' and pandoc.utils.stringify(b):match('%S') then
      local title_latex = pandoc.write(pandoc.Pandoc({ b }), 'latex')
      out:insert(pandoc.RawBlock('latex', table.concat({
        '\\begin{center}',
        '{\\Huge\\bfseries ' .. title_latex:gsub('%s+', ' ') .. '}',
        '\\end{center}',
        '\\vspace{1em}',
      }, '\n')))
      styled_title = true
    else
      if not inserted_break and b.t == 'Header' then
        out:insert(pandoc.RawBlock('latex', '\\clearpage'))
        inserted_break = true
      end
      out:insert(b)
    end
  end

  local cover = pandoc.RawBlock('latex', table.concat({
    '\\begin{center}',
    '\\includegraphics[width=0.5\\textwidth]{' .. COVER_IMAGE .. '}',
    '\\end{center}',
    '\\vspace{1em}',
  }, '\n'))
  out:insert(1, cover)

  doc.blocks = out
  return doc
end
