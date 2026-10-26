import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drawio_kit import Diagram

# Test creating the new ODD flow diagram
W = 1200
H = 1750
font_px = 14
D = Diagram(W, H, font_px=font_px)

# 1. Preparation & Inputs
D.node("p_prep", "lane", 60, 40, 1080, 140, "1. PREPARATION & INPUTS", bold_first=True, palette="lane")
D.node("inp_dtm", "data", 100, 90, 310, 70, "Digital Terrain Model\nDTM · 90 m resolution raster")
D.node("inp_canopy", "data", 445, 90, 310, 70, "Canopy Ceiling Raster\nc · 90 m resolution product (Lang et al.)")
D.node("inp_micro", "grey", 790, 90, 310, 70, "Downscaled Microclimate\nHeight-tier lookup tables (ERA5 + microclimf)")

# 2. Model Initialisation
D.node("p_init", "lane", 60, 240, 1080, 190, "2. MODEL INITIALISATION", bold_first=True, palette="lane")
D.node("init_trees", "process", 100, 290, 480, 60, "Stochastic Stand Placement\n2D Poisson process from empirical stem density")
D.node("init_cap", "process", 620, 290, 480, 60, "Bark Surface Capacity (K_xyz)\nTrunk & crown bark allometry per voxel")
D.node("init_founders", "process", 100, 360, 480, 55, "Seed Founder Individuals\nn_founders = 30 ind./taxon unit")
D.node("init_burnin", "term", 620, 360, 480, 55, "5-Generation Burn-in Period\nDiscard initial transient dynamics", bold_first=True)

# 3. Annual Simulation Loop
D.group(60, 480, 1080, 980, "3. ANNUAL SIMULATION LOOP (Years t = 1 to 50)")

# Stage 1: Reproduction & 3D Wind Dispersal
D.node("st1_box", "lane", 100, 515, 1000, 175, "STAGE 1: REPRODUCTION & 3D WIND DISPERSAL", palette="loopo")
D.node("st1_flow", "process", 140, 555, 440, 60, "Adult Flowering Eligibility\nBernoulli(p_flower(z_A)) onset gate")
D.node("st1_seeds", "process", 620, 555, 440, 60, "Seed Production\nPoisson(N · F(z)) or floor + Bernoulli residual")
D.node("st1_wind", "process", 140, 620, 440, 55, "3D Wind Dispersal Projection\nVon Mises angle + turbulent transport (Ut)")
D.node("st1_trunc", "process", 620, 620, 440, 55, "Coordinate Truncation & Boundary Check\nArray padding check (d_max, d_z)")

# Stage 2: Seedling Establishment & Space Saturation
D.node("st2_box", "lane", 100, 725, 1000, 175, "STAGE 2: SEEDLING ESTABLISHMENT & SPACE SATURATION", palette="obs")
D.node("st2_niche", "process", 140, 765, 440, 60, "Multi-Axis Climate Niche Match\nm_sp in [0, 1] suitability at landing voxel")
D.node("st2_myco", "process", 620, 765, 440, 60, "Mycorrhizal Germination Filter\nBinomial(seeds, p_germ) stochastic filter")
D.node("st2_space", "process", 140, 830, 440, 55, "Local Capacity Saturation Check\nFree space verification: C_free >= w_S")
D.node("st2_clamp", "process", 620, 830, 440, 55, "Recruit Density Clamping\nSuccessful recruits clamped to free K_xyz")

# Stage 3: Monthly Survival, Growth & Transition Loop
D.node("st3_box", "lane", 100, 935, 1000, 270, "STAGE 3: MONTHLY SURVIVAL, GROWTH & TRANSITION LOOP (m = 1 … 12)", palette="loop")
D.group(130, 970, 940, 220, "For developmental stage in {Seedling, Juvenile, Adult}:")
D.node("st3_surv", "process", 160, 1005, 880, 55, "a) Survival: Binomial(N, s_month(z, e))\nInverse-logit with beta_0 offset downscaled")
D.node("st3_trans", "process", 160, 1065, 880, 55, "b) Stage Promotion: Binomial(N_surv, q_month(e))\nTransition probability governed by psi_0 offset")
D.node("st3_growth", "process", 160, 1125, 880, 55, "c) Vegetative Size Growth: z + Delta z_m\nAdditive variance scaled pseudobulb increment")

# Stage 4: Annual Boundary Synchronization
D.node("st4_box", "lane", 100, 1240, 1000, 195, "STAGE 4: ANNUAL BOUNDARY SYNCHRONIZATION", palette="state")
D.node("st4_blend", "process", 140, 1280, 920, 45, "Seedling Cohort Size Blending: z_S(t+1) count-weighted average")
D.node("st4_cost", "process", 140, 1330, 920, 45, "Adult Reproductive Cost: 50% vegetative growth reduction on fruiting adults")
D.node("st4_update", "state", 140, 1380, 920, 45, "Update Carry-over Population Matrices for Year t + 1 (S, J, A per voxel)")

# 4. Outputs & Evaluation
D.node("p_out", "lane", 60, 1510, 1080, 190, "4. OUTPUTS & EVALUATION", bold_first=True, palette="term")
D.node("out_abund", "process", 100, 1560, 480, 60, "50-Year Abundance Trajectories\nPopulation dynamics & persistence regimes")
D.node("out_vert", "process", 620, 1560, 480, 60, "3D Realised Vertical Distribution\nCanopy profiles by life stage & taxon")
D.node("out_comp", "process", 100, 1625, 480, 60, "Interspecific Competition Displacements\nJoint vs. isolated height shifts (Kolmogorov-Smirnov D)")
D.node("out_val", "process", 620, 1625, 480, 60, "Held-Out Model Validation\nComparison against independent field records")

# Edges connecting major blocks
D.edge("inp_canopy", "init_trees", src_side="bottom", dst_side="top")
D.edge("init_founders", "st1_flow", src_side="bottom", dst_side="top", via=[(340, 450), (360, 450)])
D.edge("st1_trunc", "st2_niche", src_side="bottom", dst_side="top", via=[(840, 700), (360, 700)])
D.edge("st2_clamp", "st3_surv", src_side="bottom", dst_side="top", via=[(840, 915), (600, 915)])
D.edge("st3_growth", "st4_blend", src_side="bottom", dst_side="top", via=[(600, 1205), (600, 1270)])
D.edge("st4_update", "out_abund", src_side="bottom", dst_side="top", via=[(600, 1445), (340, 1445)], label="Year t = 50 reached", label_pos=0.5)

print("Writing files...")
D.write_drawio("/tmp/test_fig2.drawio")
D.render_png("/tmp/test_fig2.png", dpi=250)
print("Success! File size:", os.path.getsize("/tmp/test_fig2.png"))
