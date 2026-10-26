"""drawio_kit.py -- describe a flowchart once, write it both as an editable
draw.io file (.drawio, open in diagrams.net) and as a PNG rendered with
matplotlib in the same visual style.

Style follows the original hand-drawn SERES draw.io sketch (2026-09):
blue-grey process boxes, pink terminators and lane headers, yellow/orange
data stores and inputs, green loop markers and parameter inputs, purple
population state. Coordinates are draw.io pixels (origin top-left).
"""
import html, textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Ellipse, Rectangle, FancyBboxPatch, FancyArrowPatch

PAL = {
    "process":  ("#E3E9F1", "#4F6B8C"),
    "term":     ("#F8D4DA", "#B85464"),
    "data":     ("#FFF2CC", "#C9A43A"),
    "obs":      ("#FFE6CC", "#D08A1E"),
    "grey":     ("#F4F4F4", "#666666"),
    "loop":     ("#D5E8D4", "#6E9A5C"),
    "loopo":    ("#FFE6CC", "#D08A1E"),
    "param":    ("#E6F0E1", "#6E9A5C"),
    "state":    ("#E8DDF0", "#8E6AA8"),
    "decision": ("#E3E9F1", "#4F6B8C"),
    "lane":     ("#F8D4DA", "#333333"),
    "group":    ("#FBF9FE", "#8E6AA8"),
    "note":     ("#FFFFFF", "#999999"),
}
DRAWIO_SHAPE = {
    "process": "rounded=0;", "term": "ellipse;", "decision": "rhombus;",
    "data": "shape=datastore;", "grey": "shape=datastore;", "obs": "shape=parallelogram;perimeter=parallelogramPerimeter;fixedSize=1;size=12;",
    "param": "shape=parallelogram;perimeter=parallelogramPerimeter;fixedSize=1;size=12;",
    "loop": "shape=loopLimit;size=14;", "loopend": "shape=loopLimit;size=14;direction=west;",
    "loopo": "shape=loopLimit;size=14;", "state": "rounded=1;arcSize=20;", "pill": "rounded=1;arcSize=50;", "note": "rounded=1;arcSize=6;",
}


class Diagram:
    def __init__(self, width, height, font_px=11):
        self.W, self.H, self.font = width, height, font_px
        self.nodes, self.edges, self.lanes, self.groups = {}, [], [], []

    # -- building ---------------------------------------------------------
    def node(self, nid, kind, x, y, w, h, text, bold_first=True, font=None, palette=None):
        self.nodes[nid] = dict(kind=kind, x=x, y=y, w=w, h=h, text=text, bold=bold_first,
                               font=font or self.font, pal=palette or kind)
        return nid

    def lane(self, x, y, w, h, title, header=22):
        self.lanes.append(dict(x=x, y=y, w=w, h=h, title=title, header=header))

    def group(self, x, y, w, h, title, border_color=None, fill_color=None, title_pos="bottom_right", dashed=True, lw=1.2):
        self.groups.append(dict(x=x, y=y, w=w, h=h, title=title, border=border_color, fill=fill_color,
                                pos=title_pos, dashed=dashed, lw=lw))

    def edge(self, a, b, via=(), dashed=False, label=None, src_side=None, dst_side=None, label_pos=0.5,
             color=None, lw=1.0, orthogonal=True, label_color="#333333"):
        self.edges.append(dict(a=a, b=b, via=list(via), dashed=dashed, label=label,
                               src=src_side, dst=dst_side, lpos=label_pos,
                               color=color, lw=lw, orthogonal=orthogonal, lcolor=label_color))

    # -- geometry ---------------------------------------------------------
    def port(self, nid, side):
        n = self.nodes[nid]; x, y, w, h = n["x"], n["y"], n["w"], n["h"]
        return {"top": (x + w / 2, y), "bottom": (x + w / 2, y + h),
                "left": (x, y + h / 2), "right": (x + w, y + h / 2)}[side]

    def _auto_sides(self, e):
        a, b = self.nodes[e["a"]], self.nodes[e["b"]]
        ref_a = e["via"][0] if e["via"] else (b["x"] + b["w"] / 2, b["y"] + b["h"] / 2)
        ref_b = e["via"][-1] if e["via"] else (a["x"] + a["w"] / 2, a["y"] + a["h"] / 2)
        def side(n, ref):
            cx, cy = n["x"] + n["w"] / 2, n["y"] + n["h"] / 2
            dx, dy = ref[0] - cx, ref[1] - cy
            if abs(dy) * n["w"] >= abs(dx) * n["h"]:
                return "bottom" if dy > 0 else "top"
            return "right" if dx > 0 else "left"
        return e["src"] or side(a, ref_a), e["dst"] or side(b, ref_b)

    def path(self, e):
        sa, sb = self._auto_sides(e)
        p0, p1 = self.port(e["a"], sa), self.port(e["b"], sb)
        pts = [p0] + e["via"] + [p1]
        if not e.get("orthogonal", True):
            return pts, sa, sb
        # orthogonalise straight runs that are almost aligned
        out = [pts[0]]
        for q in pts[1:]:
            px, py = out[-1]
            if abs(q[0] - px) > 1 and abs(q[1] - py) > 1 and not e["via"]:
                if sa in ("top", "bottom"):
                    mid = (py + q[1]) / 2
                    out += [(px, mid), (q[0], mid)]
                else:
                    mid = (px + q[0]) / 2
                    out += [(mid, py), (mid, q[1])]
            out.append(q)
        return out, sa, sb

    # -- PNG --------------------------------------------------------------
    def _text(self, ax, n, cx, cy, w, color="#1A1A1A"):
        fs_pt = n["font"] * 0.72
        lines = n["text"].split("\n")
        chars = max(8, int((w - 10) / (0.53 * n["font"])))
        wrapped = []
        for i, ln in enumerate(lines):
            for j, piece in enumerate(textwrap.wrap(ln, chars) or [""]):
                wrapped.append((piece, i == 0 and n["bold"] and len(lines) > 1, i > 0 and n["bold"] and len(lines) > 1))
        lh = n["font"] * 1.22
        y0 = cy - (len(wrapped) - 1) * lh / 2
        for k, (t, b, sub) in enumerate(wrapped):
            ax.text(cx, y0 + k * lh, t, ha="center", va="center", fontsize=fs_pt * (0.92 if sub else 1.0),
                    fontweight="bold" if b else "normal", color="#444444" if sub else color, family="DejaVu Sans")
        return len(wrapped) * lh, max(len(t) for t, _, _ in wrapped) * 0.53 * n["font"]

    def render_png(self, path, dpi=220):
        fig = plt.figure(figsize=(self.W / 100, self.H / 100))
        ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, self.W); ax.set_ylim(self.H, 0); ax.axis("off")
        for L in self.lanes:
            f, s = PAL["lane"]
            ax.add_patch(Rectangle((L["x"], L["y"]), L["w"], L["h"], fill=False, ec=s, lw=1.0))
            ax.add_patch(Rectangle((L["x"], L["y"]), L["header"], L["h"], fc=f, ec=s, lw=1.0))
            ax.text(L["x"] + L["header"] / 2, L["y"] + L["h"] / 2, L["title"], rotation=90, ha="center", va="center",
                    fontsize=self.font * 0.72, fontweight="bold", family="DejaVu Sans")
        for G in self.groups:
            f = G.get("fill") or PAL["group"][0]
            s = G.get("border") or PAL["group"][1]
            lw = G.get("lw", 1.2)
            ls = (0, (6, 4)) if G.get("dashed", True) else "-"
            ax.add_patch(FancyBboxPatch((G["x"], G["y"]), G["w"], G["h"], boxstyle="round,pad=0,rounding_size=6",
                                        fc=f, ec=s, lw=lw, ls=ls))
            if G.get("pos") == "top_left":
                ax.text(G["x"] + 12, G["y"] + 14, G["title"], ha="left", va="top",
                        fontsize=self.font * 0.72, color=s, fontweight="bold", family="DejaVu Sans")
            else:
                ax.text(G["x"] + G["w"] - 8, G["y"] + G["h"] - 8, G["title"], ha="right", va="bottom",
                        fontsize=self.font * 0.72, color=s, fontweight="bold", family="DejaVu Sans")
        for e in self.edges:
            pts, _, _ = self.path(e)
            xs, ys = zip(*pts)
            col = e.get("color") or "#333333"
            lw = e.get("lw", 0.9)
            ls = (0, (4, 3)) if e["dashed"] else "-"
            ax.plot(xs[:-1] + (xs[-1],), ys[:-1] + (ys[-1],), color=col, lw=lw, ls=ls, solid_capstyle="butt")
            ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle="-|>", mutation_scale=9, color=col,
                                         lw=lw, shrinkA=0, shrinkB=0, linestyle="-"))
            if e["label"]:
                # place on the longest segment
                segs = list(zip(pts[:-1], pts[1:]))
                (p, q) = max(segs, key=lambda s: abs(s[1][0] - s[0][0]) + abs(s[1][1] - s[0][1]))
                lx, ly = p[0] + (q[0] - p[0]) * e["lpos"], p[1] + (q[1] - p[1]) * e["lpos"]
                ax.text(lx, ly, e["label"], ha="center", va="center", fontsize=self.font * 0.66,
                        color=e.get("lcolor", "#333333"), family="DejaVu Sans",
                        bbox=dict(fc="white", ec="none", pad=1.0))
        for nid, n in self.nodes.items():
            k = n["kind"]; f, s = PAL[n["pal"]] if n["pal"] in PAL else PAL["process"]
            x, y, w, h = n["x"], n["y"], n["w"], n["h"]
            if k == "term":
                ax.add_patch(Ellipse((x + w / 2, y + h / 2), w, h, fc=f, ec=s, lw=1.1))
            elif k == "decision":
                ax.add_patch(Polygon([(x + w / 2, y), (x + w, y + h / 2), (x + w / 2, y + h), (x, y + h / 2)], fc=f, ec=s, lw=1.1))
            elif k in ("data", "grey"):
                r = min(10, h * 0.18)
                ax.add_patch(Rectangle((x, y + r), w, h - 2 * r, fc=f, ec="none"))
                ax.add_patch(Ellipse((x + w / 2, y + h - r), w, 2 * r, fc=f, ec=s, lw=1.1))
                ax.add_patch(Rectangle((x, y + r), w, h - 2 * r, fc=f, ec="none"))
                ax.plot([x, x], [y + r, y + h - r], color=s, lw=1.1); ax.plot([x + w, x + w], [y + r, y + h - r], color=s, lw=1.1)
                for off in (0, 4):
                    ax.add_patch(Ellipse((x + w / 2, y + r + off), w, 2 * r, fc=f if off == 0 else "none", ec=s, lw=1.0))
            elif k in ("obs", "param"):
                sk = 12
                ax.add_patch(Polygon([(x + sk, y), (x + w, y), (x + w - sk, y + h), (x, y + h)], fc=f, ec=s, lw=1.1))
            elif k in ("loop", "loopo"):
                c = 12
                ax.add_patch(Polygon([(x + c, y), (x + w - c, y), (x + w, y + c), (x + w, y + h), (x, y + h), (x, y + c)], fc=f, ec=s, lw=1.1))
            elif k == "loopend":
                c = 12
                ax.add_patch(Polygon([(x, y), (x + w, y), (x + w, y + h - c), (x + w - c, y + h), (x + c, y + h), (x, y + h - c)], fc=f, ec=s, lw=1.1))
            elif k in ("state", "note", "pill"):
                r = h / 2 if k == "pill" else (14 if k == "state" else 5)
                ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=%d" % r,
                                            fc=f, ec=s, lw=1.1))
            else:
                ax.add_patch(Rectangle((x, y), w, h, fc=f, ec=s, lw=1.1))
            inner = w * (0.62 if k == "decision" else 0.9 if k in ("obs", "param") else 1.0)
            text_h, _ = self._text(ax, n, x + w / 2, y + h / 2 + (3 if k in ("data", "grey") else 0), inner)
            # usable height: diamonds and cylinders lose some of their box to the shape
            usable = h * (0.62 if k == "decision" else 1.0) - (2 * min(10, h * 0.18) if k in ("data", "grey") else 4)
            if text_h > usable:
                print(f"WARNING: text overflows node '{nid}' ({text_h:.0f} px of text in {usable:.0f} px)")
        # anything drawn outside the canvas would be cut off in the PNG
        for nid, n in self.nodes.items():
            if n["x"] < 0 or n["y"] < 0 or n["x"] + n["w"] > self.W or n["y"] + n["h"] > self.H:
                print(f"WARNING: node '{nid}' extends beyond the {self.W} x {self.H} canvas")
        fig.savefig(path, dpi=dpi, facecolor="white")
        plt.close(fig)

    # -- draw.io ----------------------------------------------------------
    def write_drawio(self, path, name="Page-1"):
        cells = ['<mxCell id="0" />', '<mxCell id="1" parent="0" />']
        def lab(t, bold):
            parts = [html.escape(p) for p in t.split("\n")]
            if bold and len(parts) > 1:
                return "<b>%s</b><br><font style=&quot;font-size: %dpx;&quot; color=&quot;#444444&quot;>%s</font>" % (
                    parts[0], self.font - 1, "<br>".join(parts[1:]))
            return "<br>".join(parts)
        for i, L in enumerate(self.lanes):
            f, s = PAL["lane"]
            cells.append(f'<mxCell id="lane{i}" value="{html.escape(L["title"])}" style="swimlane;horizontal=0;whiteSpace=wrap;html=1;'
                         f'fillColor={f};strokeColor={s};startSize={L["header"]};fontStyle=1;fontSize={self.font};" vertex="1" parent="1">'
                         f'<mxGeometry x="{L["x"]}" y="{L["y"]}" width="{L["w"]}" height="{L["h"]}" as="geometry" /></mxCell>')
        for i, G in enumerate(self.groups):
            f = G.get("fill") or PAL["group"][0]
            s = G.get("border") or PAL["group"][1]
            pos = G.get("pos", "bottom_right")
            align = "left" if pos == "top_left" else "right"
            valign = "top" if pos == "top_left" else "bottom"
            dash = "dashed=1;dashPattern=6 4;" if G.get("dashed", True) else ""
            cells.append(f'<mxCell id="grp{i}" value="{html.escape(G["title"])}" style="rounded=1;arcSize=2;whiteSpace=wrap;html=1;{dash}'
                         f'fillColor={f};strokeColor={s};fontColor={s};fontStyle=1;align={align};verticalAlign={valign};'
                         f'spacing=8;fontSize={self.font};" vertex="1" parent="1"><mxGeometry x="{G["x"]}" y="{G["y"]}" '
                         f'width="{G["w"]}" height="{G["h"]}" as="geometry" /></mxCell>')
        for nid, n in self.nodes.items():
            f, s = PAL[n["pal"]] if n["pal"] in PAL else PAL["process"]
            style = DRAWIO_SHAPE.get(n["kind"], "rounded=0;") + f"whiteSpace=wrap;html=1;fillColor={f};strokeColor={s};fontSize={n['font']};fontFamily=Helvetica;"
            cells.append(f'<mxCell id="{nid}" value="{html.escape(lab(n["text"], n["bold"]))}" style="{style}" vertex="1" parent="1">'
                         f'<mxGeometry x="{n["x"]}" y="{n["y"]}" width="{n["w"]}" height="{n["h"]}" as="geometry" /></mxCell>')
        for i, e in enumerate(self.edges):
            pts, sa, sb = self.path(e)
            ex = {"top": "0.5;exitY=0", "bottom": "0.5;exitY=1", "left": "0;exitY=0.5", "right": "1;exitY=0.5"}[sa]
            en = {"top": "0.5;entryY=0", "bottom": "0.5;entryY=1", "left": "0;entryY=0.5", "right": "1;entryY=0.5"}[sb]
            col = e.get("color") or "#333333"
            lw = e.get("lw", 1.0)
            edge_style = "edgeStyle=orthogonalEdgeStyle;" if e.get("orthogonal", True) else ""
            style = (f"{edge_style}rounded=0;html=1;endArrow=block;endFill=1;strokeColor={col};strokeWidth={lw};"
                     f"exitX={ex};entryX={en};fontSize={self.font - 1};" + ("dashed=1;" if e["dashed"] else ""))
            way = "".join(f'<mxPoint x="{p[0]}" y="{p[1]}" />' for p in pts[1:-1])
            val = html.escape(e["label"] or "")
            cells.append(f'<mxCell id="e{i}" value="{val}" style="{style}" edge="1" parent="1" source="{e["a"]}" target="{e["b"]}">'
                         f'<mxGeometry relative="1" as="geometry">' + (f'<Array as="points">{way}</Array>' if way else "") +
                         '</mxGeometry></mxCell>')
        xml = (f'<mxfile host="drawio_kit.py"><diagram name="{html.escape(name)}" id="d1"><mxGraphModel dx="1400" dy="900" grid="1" '
               f'gridSize="10" guides="1" page="1" pageWidth="{self.W}" pageHeight="{self.H}"><root>' + "".join(cells) +
               '</root></mxGraphModel></diagram></mxfile>')
        open(path, "w", encoding="utf-8").write(xml)
