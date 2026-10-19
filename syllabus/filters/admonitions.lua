-- Converts mkdocs-material style admonitions, e.g.
--
--   !!! danger "AI policy"
--
--       Body text, possibly **markdown** and [links](...).
--
-- (a paragraph beginning with "!!!" immediately followed by an indented
-- block) into styled tcolorbox environments when writing to LaTeX/PDF.
-- Passes through unchanged for all other output formats.

local style = {
  note      = {color = "gray",   label = "Note"},
  abstract  = {color = "cyan",   label = "Summary"},
  summary   = {color = "cyan",   label = "Summary"},
  info      = {color = "cyan",   label = "Info"},
  todo      = {color = "cyan",   label = "To do"},
  tip       = {color = "green",  label = "Tip"},
  hint      = {color = "green",  label = "Tip"},
  important = {color = "green",  label = "Important"},
  success   = {color = "green",  label = "Success"},
  check     = {color = "green",  label = "Success"},
  done      = {color = "green",  label = "Success"},
  question  = {color = "orange", label = "Question"},
  help      = {color = "orange", label = "Question"},
  faq       = {color = "orange", label = "Question"},
  warning   = {color = "orange", label = "Warning"},
  caution   = {color = "orange", label = "Caution"},
  attention = {color = "orange", label = "Attention"},
  failure   = {color = "red",    label = "Failure"},
  fail      = {color = "red",    label = "Failure"},
  missing   = {color = "red",    label = "Missing"},
  danger    = {color = "red",    label = "Danger"},
  error     = {color = "red",    label = "Error"},
  bug       = {color = "red",    label = "Bug"},
  example   = {color = "purple", label = "Example"},
  quote     = {color = "gray",   label = "Quote"},
  cite      = {color = "gray",   label = "Quote"},
}

local function escape_latex(s)
  s = s:gsub('\\', '\\textbackslash{}')
  s = s:gsub('([{}$&#^_%%~])', '\\%1')
  return s
end

-- Inspects the Para's inline AST directly (rather than matching on
-- stringified text) so a quoted title such as "AI policy" is recognized
-- even though pandoc's smart-typography reader turns it into a Quoted
-- inline (curly quotes) rather than a literal Str with '"' characters.
local function admonition_kind_and_title(para)
  local c = para.content
  local idx = 1

  if not (c[idx] and c[idx].t == 'Str' and c[idx].text == '!!!') then
    return nil, nil
  end
  idx = idx + 1

  while c[idx] and c[idx].t == 'Space' do idx = idx + 1 end
  local kind_tok = c[idx]
  if not (kind_tok and kind_tok.t == 'Str') then return nil, nil end
  local kind = kind_tok.text:lower()
  idx = idx + 1

  while c[idx] and c[idx].t == 'Space' do idx = idx + 1 end
  local title = nil
  if c[idx] and c[idx].t == 'Quoted' then
    title = pandoc.utils.stringify(pandoc.Inlines(c[idx].content))
    idx = idx + 1
  end

  while c[idx] and c[idx].t == 'Space' do idx = idx + 1 end
  if c[idx] then return nil, nil end -- trailing content: not an admonition marker

  return kind, title
end

function Blocks(blocks)
  if FORMAT ~= 'latex' and FORMAT ~= 'beamer' then
    return blocks
  end

  local out = pandoc.Blocks({})
  local i = 1
  local n = #blocks

  while i <= n do
    local b = blocks[i]
    local next_b = blocks[i + 1]

    if b.t == 'Para' and next_b and next_b.t == 'CodeBlock' then
      local kind, title = admonition_kind_and_title(b)
      local s = kind and style[kind]

      if s then
        local heading = (title and title ~= '') and title or s.label
        local inner_blocks = pandoc.read(next_b.text, 'markdown').blocks
        local body_latex = pandoc.write(pandoc.Pandoc(inner_blocks), 'latex')
        local box = string.format(
          '\\begin{tcolorbox}[breakable,colback=%s!8,colframe=%s!65!black,title=%s,fonttitle=\\bfseries]\n%s\\end{tcolorbox}',
          s.color, s.color, escape_latex(heading), body_latex)
        out:insert(pandoc.RawBlock('latex', box))
        i = i + 2
        goto continue
      end
    end

    out:insert(b)
    i = i + 1
    ::continue::
  end

  return out
end
