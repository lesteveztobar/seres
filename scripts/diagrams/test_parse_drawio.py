import xml.etree.ElementTree as ET
import html, re

def parse_val(val):
    if not val:
        return ""
    # replace breaks and div closings with newline
    s = val.replace("<br>", "\n").replace("<br/>", "\n").replace("<br />", "\n")
    s = s.replace("</div>", "\n").replace("</p>", "\n")
    # remove html tags
    s = re.sub(r'<[^>]+>', '', s)
    s = html.unescape(s)
    # clean extra empty lines
    lines = [ln.strip() for ln in s.split("\n")]
    # remove leading/trailing empty lines
    while lines and not lines[0]:
        lines.pop(0)
    while lines and not lines[-1]:
        lines.pop(-1)
    return "\n".join(lines)

for f in ['report/figures/fig_causal_diagram.drawio', 'report/figures/fig2_seres_sketch.drawio', 'report/figures/fig_lookup_climate.drawio']:
    print("=== Testing", f, "===")
    tree = ET.parse(f)
    root = tree.getroot()
    model = root.find('.//mxGraphModel')
    print("Canvas:", model.get('pageWidth'), "x", model.get('pageHeight'))
    count_v = 0
    count_e = 0
    for cell in root.findall('.//mxCell'):
        if cell.get('vertex'):
            count_v += 1
            t = parse_val(cell.get('value', ''))
            if t:
                first_ln = t.split('\n')[0]
                # print(f"  Node {cell.get('id')}: {first_ln[:40]}")
        elif cell.get('edge'):
            count_e += 1
    print(f"  Total vertices: {count_v}, edges: {count_e}")
print("Done!")
