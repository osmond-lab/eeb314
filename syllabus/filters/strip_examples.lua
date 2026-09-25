-- Drops "- [Example](...)" list items (links to final_project/partI_example.md
-- etc.) from the PDF. Those links only make sense on the live website;
-- the linked .md files aren't part of the printed syllabus.

local function is_example_item(item)
  if #item ~= 1 then return false end
  local b = item[1]
  if b.t ~= 'Plain' and b.t ~= 'Para' then return false end
  if #b.content ~= 1 then return false end
  local inline = b.content[1]
  if inline.t ~= 'Link' then return false end
  return pandoc.utils.stringify(inline.content) == 'Example'
end

function BulletList(el)
  if FORMAT ~= 'latex' and FORMAT ~= 'beamer' then
    return el
  end
  local kept = {}
  for _, item in ipairs(el.content) do
    if not is_example_item(item) then
      table.insert(kept, item)
    end
  end
  el.content = kept
  return el
end
