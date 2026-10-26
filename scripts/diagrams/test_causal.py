import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drawio_kit import Diagram

# Diagram: Biological Causal Diagram (Eco-Physiological Relationships)
W = 1200
H = 1600
font_px = 14
D = Diagram(W, H, font_px=font_px)

# 1. Synoptic Macroclimate (ERA5)
D.node("macro", "process", 250, 40, 700, 75,
       "SYNOPTIC MACROCLIMATE (ERA5)\nTemperature, Relative Humidity, Radiation, Wind",
       bold_first=True, palette="grey")

# 2. 3D Voxel Canopy Microclimate
D.node("micro", "process", 200, 180, 800, 85,
       "3D VOXEL CANOPY MICROCLIMATE\nDaytime / Night Temperature (T), Relative Humidity (RH),\nSolar Radiation (S), Effective Wind Field (u)",
       bold_first=True, palette="process")

# Branch left: Substrate Canopy Geometry & Bark
D.node("geom", "data", 80, 340, 340, 105,
       "SUBSTRATE CANOPY GEOMETRY & BARK\nAllometric Trunk & Crown Bark Area\nCarrying Capacity (K_xyz) per voxel",
       bold_first=True, palette="data")

# Branch right: Niche Suitability
D.node("niche", "obs", 780, 340, 340, 105,
       "NICHE SUITABILITY\n4-Axis Multiplicative KDE Intersection\nm_sp(x) in [0, 1]",
       bold_first=True, palette="obs")

# Center: Seedling Establishment
D.node("est", "term", 400, 510, 400, 95,
       "SEEDLING ESTABLISHMENT\nSpace Clamping & Mycorrhizal Gate\n(p_germ ~ 10^-3, C_free >= w_S)",
       bold_first=True, palette="term")

# Demographic Life Cycle Group
D.group(80, 670, 1040, 470, "DEMOGRAPHIC LIFE CYCLE")

# Seedlings
D.node("seedlings", "loop", 140, 720, 260, 75,
       "SEEDLINGS\nz in [0.2, 1.0] cm\nCohort size blending z_S",
       bold_first=True, palette="loop")

# Juveniles
D.node("juveniles", "loopo", 500, 720, 260, 75,
       "JUVENILES\nz in [1.0, 7.0] cm\nIntermediate vegetative stage",
       bold_first=True, palette="loopo")

# Adults
D.node("adults", "state", 820, 870, 260, 80,
       "ADULTS\nz >= 7.0 cm\nFlowering & fruit onset",
       bold_first=True, palette="state")

# Growth & Survival boxes
D.node("surv_s", "process", 140, 860, 260, 60,
       "Voxel Survival & Growth\ns_month(z, e), Delta z_m",
       bold_first=False, palette="process")

D.node("surv_j", "process", 500, 860, 260, 60,
       "Voxel Survival & Growth\ns_month(z, e), Delta z_m",
       bold_first=False, palette="process")

D.node("surv_a", "process", 820, 1010, 260, 60,
       "Adult Survival & Maintenance\ns_month(z, e), 50% fruit cost",
       bold_first=False, palette="process")

# 3D Wind Dispersal Kernel
D.node("wind", "process", 280, 1220, 640, 95,
       "3D WIND DISPERSAL KERNEL\nTerminal Velocity (Ut = 0.23 m/s) · Turbulent Airflow (lambda = 3.23)\nBounded Dispersal Window (d_max, d_z)",
       bold_first=True, palette="loopo")

# Edges
D.edge("macro", "micro", label="Topographic & Canopy Attenuation", label_pos=0.5)
D.edge("micro", "geom", src_side="bottom", dst_side="top", via=[(350, 300), (250, 300)])
D.edge("micro", "niche", src_side="bottom", dst_side="top", via=[(850, 300), (950, 300)])
D.edge("micro", "est", src_side="bottom", dst_side="top", label="Microclimate Gating", label_pos=0.5)

D.edge("geom", "est", src_side="bottom", dst_side="left", via=[(250, 555)], label="Space Saturation (C_free >= w_S)", label_pos=0.5)
D.edge("niche", "est", src_side="bottom", dst_side="right", via=[(950, 555)], label="Environmental Gating", label_pos=0.5)

D.edge("est", "seedlings", src_side="bottom", dst_side="top", via=[(600, 640), (270, 640)], label="Recruits into S cohort", label_pos=0.6)

D.edge("seedlings", "juveniles", label="Promote (psi)", label_pos=0.5)
D.edge("seedlings", "surv_s")
D.edge("juveniles", "adults", src_side="right", dst_side="top", via=[(800, 755), (950, 755)], label="Promote (psi)", label_pos=0.5)
D.edge("juveniles", "surv_j")
D.edge("adults", "surv_a")

D.edge("adults", "wind", src_side="bottom", dst_side="top", via=[(950, 1170), (600, 1170)],
       label="Reproductive Onset & Seed Supply: p_flower(z_A), n_fruits x S", label_pos=0.5)
D.edge("wind", "est", src_side="left", dst_side="left", via=[(120, 1267), (120, 555)], dashed=True,
       label="Seed Rain to Candidate Voxels", label_pos=0.5)

print("Writing files...")
D.write_drawio("/tmp/test_causal.drawio")
D.render_png("/tmp/test_causal.png", dpi=250)
print("Success! File size:", os.path.getsize("/tmp/test_causal.png"))
