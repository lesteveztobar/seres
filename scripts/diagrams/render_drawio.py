"""render_drawio.py -- renders draw.io XML files to crisp high-resolution PNGs
using matplotlib, preserving exact node geometries, colors, typography, and edge routing.
"""
import os, sys, re, html
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Ellipse, Rectangle, FancyBboxPatch, FancyArrowPatch
import textwrap

def parse_style(s):
    if not s:
        return {}
    res = {}
    for part in s.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            res[k.strip()] = v.strip()
        else:
            res[part] = True
    return res

def clean_html_text(val):
    if not val:
        return []
    # Replace line breaks and paragraph/div endings with newline
    s = val.replace("<br>", "\n").replace("<br/>", "\n").replace("<br />", "\n")
    s = s.replace("</div>", "\n").replace("</p>", "\n")
    
    # Check for bold tag on first line or pieces
    lines = s.split("\n")
    cleaned_lines = []
    for ln in lines:
        ln_strip = ln.strip()
        if not ln_strip:
            continue
        is_bold = ("<b>" in ln or "font-weight: bold" in ln or "font-weight:bold" in ln or "fontStyle=1" in ln)
        is_italic = ("<i>" in ln or "font-style: italic" in ln)
        # extract font color if present
        color_match = re.search(r'color=["\']?#?([0-9a-fA-F]{6}|[0-9a-fA-F]{3})["\']?', ln)
        font_color = ("#" + color_match.group(1)) if color_match else None
        
        # strip tags
        txt = re.sub(r'<[^>]+>', '', ln)
        txt = html.unescape(txt).strip()
        if txt:
            cleaned_lines.append(dict(text=txt, bold=is_bold, italic=is_italic, color=font_color))
    return cleaned_lines

def render_drawio(drawio_path, png_path, dpi=300):
    tree = ET.parse(drawio_path)
    root = tree.getroot()
    model = root.find('.//mxGraphModel')
    page_w = float(model.get('pageWidth', 1100))
    page_h = float(model.get('pageHeight', 800))
    
    vertices = []
    edges = []
    id_map = {}
    
    # 1. Parse all cells
    for cell in root.findall('.//mxCell'):
        cid = cell.get('id')
        style = parse_style(cell.get('style', ''))
        val = cell.get('value', '')
        geo = cell.find('mxGeometry')
        
        if cell.get('vertex') and geo is not None:
            x = float(geo.get('x', 0))
            y = float(geo.get('y', 0))
            w = float(geo.get('width', 0))
            h = float(geo.get('height', 0))
            node_data = dict(id=cid, x=x, y=y, w=w, h=h, style=style, val=val, cell=cell)
            vertices.append(node_data)
            id_map[cid] = node_data
            
        elif cell.get('edge'):
            src = cell.get('source')
            dst = cell.get('target')
            pts = []
            src_pt = None
            dst_pt = None
            if geo is not None:
                for pt in geo.findall('.//mxPoint'):
                    pt_type = pt.get('as')
                    px = float(pt.get('x', 0))
                    py = float(pt.get('y', 0))
                    if pt_type == 'sourcePoint':
                        src_pt = (px, py)
                    elif pt_type == 'targetPoint':
                        dst_pt = (px, py)
                    else:
                        pts.append((px, py))
            edges.append(dict(id=cid, src=src, dst=dst, pts=pts, src_pt=src_pt, dst_pt=dst_pt,
                              style=style, val=val, cell=cell))
            
    # 2. Compute bounding box using true oriented bounds
    min_x = 1e9
    min_y = 1e9
    max_x = -1e9
    max_y = -1e9
    
    for v in vertices:
        rot = float(v['style'].get('rotation', 0))
        cx = v['x'] + v['w'] / 2
        cy = v['y'] + v['h'] / 2
        if rot in (90, -90, 270):
            w_eff, h_eff = v['h'], v['w']
        else:
            w_eff, h_eff = v['w'], v['h']
        min_x = min(min_x, cx - w_eff / 2)
        min_y = min(min_y, cy - h_eff / 2)
        max_x = max(max_x, cx + w_eff / 2)
        max_y = max(max_y, cy + h_eff / 2)

    for e in edges:
        for pt in e['pts']:
            min_x = min(min_x, pt[0])
            min_y = min(min_y, pt[1])
            max_x = max(max_x, pt[0])
            max_y = max(max_y, pt[1])
        
    pad = 35
    canvas_x0 = min_x - pad
    canvas_y0 = min_y - pad
    canvas_w = (max_x - min_x) + 2 * pad
    canvas_h = (max_y - min_y) + 2 * pad
    
    fig = plt.figure(figsize=(canvas_w / 100, canvas_h / 100))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(canvas_x0, canvas_x0 + canvas_w)
    ax.set_ylim(canvas_y0 + canvas_h, canvas_y0) # origin top-left
    ax.axis("off")
    
    # 3. Sort vertices: backgrounds/containers first, foreground nodes second
    # Containers have larger area or are dashed groups
    def vertex_order(v):
        area = v['w'] * v['h']
        is_container = v['style'].get('dashed') or area > 70000 or 'swimlane' in v['style']
        return (0 if is_container else 1, -area)
        
    sorted_vertices = sorted(vertices, key=vertex_order)
    
    # Render background shapes
    for v in sorted_vertices:
        st = v['style']
        area = v['w'] * v['h']
        is_container = st.get('dashed') or area > 70000 or 'swimlane' in st
        if not is_container:
            continue
        _draw_shape(ax, v)

    # 4. Render edges
    for e in edges:
        _draw_edge(ax, e, id_map)
        
    # 5. Render foreground nodes
    for v in sorted_vertices:
        st = v['style']
        area = v['w'] * v['h']
        is_container = st.get('dashed') or area > 70000 or 'swimlane' in st
        if is_container:
            continue
        _draw_shape(ax, v)
        
    # 6. Render labels for all nodes
    for v in sorted_vertices:
        _draw_labels(ax, v)

    # 7. Render edge labels
    for e in edges:
        _draw_edge_label(ax, e, id_map)

    os.makedirs(os.path.dirname(png_path), exist_ok=True)
    fig.savefig(png_path, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"Rendered: {png_path} ({canvas_w:.0f}x{canvas_h:.0f} px at {dpi} dpi)")

def _draw_shape(ax, v):
    st = v['style']
    x, y, w, h = v['x'], v['y'], v['w'], v['h']
    
    fc = st.get('fillColor', '#FFFFFF')
    sc = st.get('strokeColor', '#333333')
    if fc in ('none', 'transparent', None):
        fc = 'none'
    if sc in ('none', 'transparent', None):
        sc = 'none'
        
    lw = float(st.get('strokeWidth', 1.0))
    is_dashed = (st.get('dashed') in ('1', True) or 'dashPattern' in st)
    ls = (0, (6, 4)) if is_dashed else '-'
    
    is_text = ('text' in st or st.get('shape') == 'text')
    if is_text:
        return
        
    is_ellipse = ('ellipse' in st or st.get('shape') == 'ellipse')
    is_rhombus = ('rhombus' in st or st.get('shape') == 'rhombus')
    is_datastore = (st.get('shape') == 'datastore')
    is_rounded = (st.get('rounded') in ('1', True))
    arc_size = float(st.get('arcSize', 20 if is_rounded else 0))
    
    if is_ellipse:
        ax.add_patch(Ellipse((x + w / 2, y + h / 2), w, h, fc=fc, ec=sc, lw=lw, ls=ls))
    elif is_rhombus:
        pts = [(x + w / 2, y), (x + w, y + h / 2), (x + w / 2, y + h), (x, y + h / 2)]
        ax.add_patch(Polygon(pts, fc=fc, ec=sc, lw=lw, ls=ls))
    elif is_datastore:
        r = min(12, h * 0.18)
        # body
        ax.add_patch(Rectangle((x, y + r), w, h - 2 * r, fc=fc, ec="none"))
        # bottom rim
        ax.add_patch(Ellipse((x + w / 2, y + h - r), w, 2 * r, fc=fc, ec=sc, lw=lw))
        ax.add_patch(Rectangle((x, y + r), w, h - 2 * r, fc=fc, ec="none"))
        ax.plot([x, x], [y + r, y + h - r], color=sc, lw=lw)
        ax.plot([x + w, x + w], [y + r, y + h - r], color=sc, lw=lw)
        # top rim
        for off in (0, 3):
            ax.add_patch(Ellipse((x + w / 2, y + r + off), w, 2 * r, fc=fc if off == 0 else "none", ec=sc, lw=lw))
    elif is_rounded:
        r_size = h / 2 if arc_size >= 40 else min(14, arc_size * 0.5)
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r_size}",
                                    fc=fc, ec=sc, lw=lw, ls=ls))
    else:
        ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=sc, lw=lw, ls=ls))

def _draw_edge(ax, e, id_map):
    st = e['style']
    src_id = e['src']
    dst_id = e['dst']
    
    src = id_map.get(src_id)
    dst = id_map.get(dst_id)
    
    # exit/entry fraction
    exit_x = float(st.get('exitX', 0.5))
    exit_y = float(st.get('exitY', 0.5))
    entry_x = float(st.get('entryX', 0.5))
    entry_y = float(st.get('entryY', 0.5))
    
    if e['src_pt']:
        p0 = e['src_pt']
    elif src:
        p0 = (src['x'] + src['w'] * exit_x, src['y'] + src['h'] * exit_y)
    else:
        p0 = e['pts'][0] if e['pts'] else (0, 0)
        
    if e['dst_pt']:
        p1 = e['dst_pt']
    elif dst:
        p1 = (dst['x'] + dst['w'] * entry_x, dst['y'] + dst['h'] * entry_y)
    else:
        p1 = e['pts'][-1] if e['pts'] else (0, 0)
        
    all_pts = [p0] + e['pts'] + [p1]
    
    # If edge is orthogonal and only 2 points given with straight angle
    if st.get('edgeStyle') == 'orthogonalEdgeStyle' and not e['pts'] and src and dst:
        # compute auto mid
        dx = abs(p1[0] - p0[0])
        dy = abs(p1[1] - p0[1])
        if dx > 1 and dy > 1:
            if exit_y in (0, 1) or entry_y in (0, 1):
                mid_y = (p0[1] + p1[1]) / 2
                all_pts = [p0, (p0[0], mid_y), (p1[0], mid_y), p1]
            else:
                mid_x = (p0[0] + p1[0]) / 2
                all_pts = [p0, (mid_x, p0[1]), (mid_x, p1[1]), p1]

    xs, ys = zip(*all_pts)
    col = st.get('strokeColor', '#333333')
    if col in ('none', 'transparent', None):
        col = '#333333'
    lw = float(st.get('strokeWidth', 1.0))
    is_dashed = (st.get('dashed') in ('1', True) or 'dashPattern' in st)
    ls = (0, (4, 3)) if is_dashed else '-'
    
    ax.plot(xs, ys, color=col, lw=lw, ls=ls, solid_capstyle="butt", zorder=2)
    # arrow head
    if st.get('endArrow') != 'none' and len(all_pts) >= 2:
        ax.add_patch(FancyArrowPatch(all_pts[-2], all_pts[-1], arrowstyle="-|>",
                                     mutation_scale=10, color=col, lw=lw,
                                     shrinkA=0, shrinkB=0, linestyle="-", zorder=3))

def _draw_labels(ax, v):
    st = v['style']
    x, y, w, h = v['x'], v['y'], v['w'], v['h']
    rot = float(st.get('rotation', 0))
    
    cleaned = clean_html_text(v['val'])
    if not cleaned:
        return
        
    base_fs = float(st.get('fontSize', 13))
    align = st.get('align', 'center')
    valign = st.get('verticalAlign', 'middle')
    font_family = "DejaVu Sans"
    
    cx = x + w / 2
    cy = y + h / 2
    
    # Handle rotated text (e.g. swimlane vertical headers)
    if rot in (90, -90, 270):
        # matplotlib rotation: 90 is counter-clockwise (bottom to top), -90/270 is top to bottom
        plot_rot = 90 if rot in (-90, 270) else -90
        # draw at center of node box
        full_text = "\n".join(item['text'] for item in cleaned)
        col = st.get('fontColor') or cleaned[0].get('color') or '#333333'
        ax.text(cx, cy, full_text, rotation=plot_rot, ha="center", va="center",
                fontsize=base_fs * 0.72, fontweight="bold", fontstyle="italic",
                color=col, family=font_family, zorder=5)
        return

    # Normal text
    lh = base_fs * 1.25
    total_h = len(cleaned) * lh
    
    if valign == 'top':
        start_y = y + 14 + lh / 2
    elif valign == 'bottom':
        start_y = y + h - 14 - total_h + lh / 2
    else:
        start_y = cy - (total_h - lh) / 2
        
    ha = align if align in ('left', 'right') else 'center'
    text_x = (x + 14) if ha == 'left' else ((x + w - 14) if ha == 'right' else cx)
    
    # Wrap lines that exceed the usable box width
    chars_per_line = max(10, int((w - 16) / (0.52 * base_fs)))
    wrapped_items = []
    for item in cleaned:
        txt = item['text']
        if len(txt) > chars_per_line:
            pieces = textwrap.wrap(txt, width=chars_per_line) or [txt]
            for p in pieces:
                wrapped_items.append(dict(text=p, bold=item['bold'], italic=item['italic'], color=item['color']))
        else:
            wrapped_items.append(item)

    lh = base_fs * 1.25
    total_h = len(wrapped_items) * lh
    
    if valign == 'top':
        start_y = y + 14 + lh / 2
    elif valign == 'bottom':
        start_y = y + h - 14 - total_h + lh / 2
    else:
        start_y = cy - (total_h - lh) / 2
    
    for i, item in enumerate(wrapped_items):
        cur_y = start_y + i * lh
        col = item.get('color') or st.get('fontColor') or ('#1A1A1A' if (item['bold'] or i == 0) else '#444444')
        fs = (base_fs if (item['bold'] or i == 0) else base_fs * 0.92) * 0.72
        fw = 'bold' if item['bold'] else 'normal'
        fst = 'italic' if item['italic'] else 'normal'
        ax.text(text_x, cur_y, item['text'], ha=ha, va="center",
                fontsize=fs, fontweight=fw, fontstyle=fst,
                color=col, family=font_family, zorder=5)

def _draw_edge_label(ax, e, id_map):
    st = e['style']
    val = e['val']
    if not val:
        return
    cleaned = clean_html_text(val)
    if not cleaned:
        return
    label_txt = "\n".join(item['text'] for item in cleaned)
    
    # get coordinates
    src_id = e['src']
    dst_id = e['dst']
    src = id_map.get(src_id)
    dst = id_map.get(dst_id)
    exit_x = float(st.get('exitX', 0.5))
    exit_y = float(st.get('exitY', 0.5))
    entry_x = float(st.get('entryX', 0.5))
    entry_y = float(st.get('entryY', 0.5))
    
    p0 = e['src_pt'] or ((src['x'] + src['w'] * exit_x, src['y'] + src['h'] * exit_y) if src else (0, 0))
    p1 = e['dst_pt'] or ((dst['x'] + dst['w'] * entry_x, dst['y'] + dst['h'] * entry_y) if dst else (0, 0))
    all_pts = [p0] + e['pts'] + [p1]
    
    # place on the longest segment
    segs = list(zip(all_pts[:-1], all_pts[1:]))
    (p, q) = max(segs, key=lambda s: abs(s[1][0] - s[0][0]) + abs(s[1][1] - s[0][1]))
    lx = (p[0] + q[0]) / 2
    ly = (p[1] + q[1]) / 2
    
    fs = float(st.get('fontSize', 12)) * 0.72
    ax.text(lx, ly, label_txt, ha="center", va="center", fontsize=fs,
            family="DejaVu Sans", color="#333333",
            bbox=dict(fc="white", ec="none", pad=1.2), zorder=6)

def main():
    files = [
        ("report/figures/fig_causal_diagram.drawio", "report/figures/fig_causal_diagram.png"),
        ("report/figures/fig2_seres_sketch.drawio", "report/figures/fig2_seres_sketch.png"),
        ("report/figures/fig_lookup_climate.drawio", "report/figures/fig_lookup_climate.png"),
    ]
    for src, dst in files:
        if os.path.exists(src):
            print(f"Processing: {src} -> {dst}")
            render_drawio(src, dst, dpi=300)
            # Also mirror to manuscript for paper/figures/
            manuscript_dst = os.path.join("manuscript for paper", "figures", os.path.basename(dst))
            manuscript_src = os.path.join("manuscript for paper", "figures", os.path.basename(src))
            # copy drawio and png to manuscript
            with open(src, "r", encoding="utf-8") as f_in, open(manuscript_src, "w", encoding="utf-8") as f_out:
                f_out.write(f_in.read())
            render_drawio(src, manuscript_dst, dpi=300)
        else:
            print(f"File not found: {src}")

if __name__ == "__main__":
    main()
