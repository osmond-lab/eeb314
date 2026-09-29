-- Drops the bulleted list of textbook links (e-copies, physical copies,
-- buy your own) from the printed PDF -- those links only make sense on
-- the live website. The paragraph above (author/title) is kept.

local in_textbook = false

function Header(el)
  in_textbook = pandoc.utils.stringify(el.content) == 'Textbook'
  return el
end

function BulletList(el)
  if FORMAT ~= 'latex' and FORMAT ~= 'beamer' then
    return el
  end
  if not in_textbook then
    return el
  end
  return {}
end
