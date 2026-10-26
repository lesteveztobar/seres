"""make_diagrams.py -- builds the flowcharts and conceptual diagrams of the thesis
from drawio_kit.py, each as an editable .drawio file and a crisp high-res PNG:
  fig2_seres_sketch.{drawio,png}      (Program Execution Flow Diagram / ODD Scheduling Sequence)
  fig_causal_diagram.{drawio,png}     (Biological Causal Diagram / Eco-Physiological Network)
  fig_lookup_climate.{drawio,png}     (Climate Lookup by Height, Appendix)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drawio_kit import Diagram

OUT = os.environ.get("DIAGRAM_OUT") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "report", "figures")


def biological_causal_diagram():
    """Biological Causal Diagram (Eco-Physiological Network)
    Faithful to user prototype sketch #1 (media_...0ba6da52.png).
    Zero programming jargon -- purely biological and eco-physiological concepts.
    """
    W = 1120
    H = 760
    font_px = 13.5
    D = Diagram(W, H, font_px=font_px)

    # Top: System Outcome
    D.node("occ", "pill", 300, 25, 520, 58,
           "Species Occurrence & Abundance\nCanopy distribution, population presence, and local density",
           bold_first=True, palette="term")

    # Upper demographic endpoints / processes
    D.node("poaching", "term", 45, 125, 185, 54,
           "Poaching\nAnthropogenic harvesting",
           bold_first=True, palette="grey")

    D.node("dispersal", "term", 280, 125, 195, 54,
           "Dispersal\n3D canopy seed transport",
           bold_first=True, palette="process")

    D.node("establishment", "term", 575, 125, 215, 54,
           "Establishment\nGermination & anchor on bark",
           bold_first=True, palette="process")

    D.node("reproduction", "term", 870, 205, 195, 54,
           "Reproduction\nFruit & viable seed set",
           bold_first=True, palette="term")

    # Middle layer: Vital rates & demographics
    D.node("survival", "term", 385, 215, 205, 54,
           "Survival\nPersistence across life stages",
           bold_first=True, palette="process")

    D.node("growth", "term", 270, 325, 195, 54,
           "Growth\nBiomass & vegetative size",
           bold_first=True, palette="process")

    D.node("pollination", "term", 645, 285, 195, 54,
           "Pollination\nPollen transfer & fertilization",
           bold_first=True, palette="process")

    D.node("flowering", "term", 875, 335, 195, 54,
           "Flowering\nInflorescence onset & phenology",
           bold_first=True, palette="process")

    D.node("competition", "term", 455, 425, 185, 54,
           "Competition\nBark space crowding",
           bold_first=True, palette="process")

    D.node("mortality", "term", 665, 425, 185, 54,
           "Mortality\nDesiccation & senescence",
           bold_first=True, palette="process")

    # Foundational Environmental Drivers (Soft green)
    D.node("connectivity", "term", 45, 285, 195, 54,
           "Canopy Connectivity\nCanopy bridges & wind corridors",
           bold_first=True, palette="param")

    D.node("nutrition", "term", 155, 510, 215, 54,
           "Substrate Nutrition\nMoss mats & moisture retention",
           bold_first=True, palette="param")

    D.node("habitat", "term", 445, 630, 245, 54,
           "Host Tree Habitat\nTree architecture & bark capacity",
           bold_first=True, palette="param")

    D.node("microclimate", "term", 745, 520, 260, 54,
           "Canopy Microenvironment\nStratified T, RH, and solar buffering",
           bold_first=True, palette="param")

    # Edges (Direct straight arrows matching sketch #1)
    # 1. Poaching impacts
    D.edge("poaching", "occ", orthogonal=False, lw=1.1)
    D.edge("poaching", "survival", orthogonal=False, lw=1.1)

    # 2. Recruitment & occurrence
    D.edge("establishment", "occ", orthogonal=False, lw=1.2)
    D.edge("dispersal", "establishment", orthogonal=False, lw=1.2)
    D.edge("establishment", "growth", orthogonal=False, lw=1.1, dashed=True)

    # 3. Reproductive cascade
    D.edge("reproduction", "dispersal", orthogonal=False, lw=1.2)
    D.edge("pollination", "reproduction", orthogonal=False, lw=1.2)
    D.edge("flowering", "pollination", orthogonal=False, lw=1.2)
    D.edge("survival", "reproduction", orthogonal=False, lw=1.1)

    # 4. Growth & survival interactions
    D.edge("growth", "survival", orthogonal=False, lw=1.2)
    D.edge("growth", "flowering", orthogonal=False, lw=1.1)
    D.edge("competition", "survival", orthogonal=False, lw=1.1)
    D.edge("mortality", "survival", orthogonal=False, lw=1.1)

    # 5. Connectivity drivers
    D.edge("connectivity", "dispersal", orthogonal=False, lw=1.1, color="#2E7D32")
    D.edge("connectivity", "pollination", orthogonal=False, lw=1.1, color="#2E7D32")

    # 6. Environmental cascades
    D.edge("nutrition", "growth", orthogonal=False, lw=1.1, color="#2E7D32")
    D.edge("habitat", "nutrition", orthogonal=False, lw=1.1, color="#2E7D32")
    D.edge("habitat", "microclimate", orthogonal=False, lw=1.1, color="#2E7D32")

    # 7. Microclimate influences
    D.edge("microclimate", "mortality", orthogonal=False, lw=1.1, color="#2E7D32")
    D.edge("microclimate", "competition", orthogonal=False, lw=1.1, color="#2E7D32")
    D.edge("microclimate", "flowering", orthogonal=False, lw=1.1, color="#2E7D32")
    D.edge("microclimate", "establishment", orthogonal=False, lw=1.1, color="#2E7D32")
    D.edge("microclimate", "survival", orthogonal=False, lw=1.1, color="#2E7D32")

    return D


def odd_execution_flow():
    """Program Execution Flow Diagram (ODD Scheduling Sequence)
    Vertical A4 format matching user prototype sketch #2 (media_...5862aefb.png).
    Clearly organizes the simulation hierarchy:
      1. Preparation & Inputs
      2. Model Initialisation
      3. Annual Generational Loop enclosing Monthly Demographics Loop
      4. Outputs & Virtual Ecologist Observation Module
    """
    W = 1120
    H = 1720
    font_px = 13.5
    D = Diagram(W, H, font_px=font_px)

    # 1. Preparation & Inputs
    D.node("p_prep", "lane", 50, 25, 1020, 135, "1. PREPARATION & INPUT DATA", bold_first=True, palette="lane")
    D.node("inp_dtm", "data", 80, 68, 295, 76,
           "Digital Terrain Model\nDTM · 90 m raster topography")
    D.node("inp_canopy", "data", 410, 68, 295, 76,
           "Canopy Ceiling Raster\nc · 90 m canopy top (Lang et al.)")
    D.node("inp_micro", "grey", 740, 68, 300, 76,
           "Downscaled Microclimate\nHeight-tier lookup tables (ERA5 + microclimf)")

    # 2. Model Initialisation
    D.node("p_init", "lane", 50, 180, 1020, 175, "2. MODEL INITIALISATION", bold_first=True, palette="lane")
    D.node("init_trees", "process", 80, 225, 465, 55,
           "Stochastic Host Tree Placement\n2D Poisson process from empirical stem density")
    D.node("init_cap", "process", 575, 225, 465, 55,
           "Bark Surface Capacity Allocation (K_xyz)\nTrunk & crown bark allometry per 3D voxel")
    D.node("init_founders", "process", 80, 290, 465, 50,
           "Seed Founder Individuals\nn_founders = 30 ind./taxon unit seeded into domain")
    D.node("init_burnin", "term", 575, 290, 465, 50,
           "5-Generation Burn-in Period\nDiscard initial transient dynamics before tracking",
           bold_first=True)

    # 3. Annual Simulation Loop (Outer Purple Container)
    D.group(50, 375, 1020, 1020,
            "3. ANNUAL SIMULATION LOOP (Years t = 1 … 50)  [1 generation / timestep]",
            border_color="#8E6AA8", fill_color="#FBF9FE", title_pos="top_left", lw=1.5)

    # Environmental Driver Nodes on Left
    D.node("env_hab", "state", 80, 440, 220, 52,
           "Host Habitat Structure\nTree crown allometry & K_xyz",
           bold_first=True, palette="param")
    D.node("env_micro", "state", 80, 540, 220, 58,
           "Canopy Microclimate\nT, RH, radiation height tiers",
           bold_first=True, palette="param")
    D.node("env_nutr", "state", 80, 755, 220, 52,
           "Substrate Nutrition\nMoisture & moss cover",
           bold_first=True, palette="param")
    D.node("env_conn", "state", 80, 935, 220, 52,
           "Canopy Connectivity\nAerodynamic wind corridors",
           bold_first=True, palette="param")

    # Inner Orange Container: Monthly Demographics Loop
    D.group(335, 420, 705, 295,
            "STAGE 3: MONTHLY DEMOGRAPHICS LOOP (m = 1 … 12)  [1 month / timestep]",
            border_color="#E67E22", fill_color="#FFFDF9", title_pos="top_left", lw=1.4)

    D.node("st3_growth", "process", 360, 460, 430, 55,
           "Vegetative Size Growth (z + Δz_m)\nAdditive variance scaled pseudobulb increment",
           bold_first=True, palette="process")

    D.node("st3_surv", "process", 360, 530, 430, 80,
           "Survival & Stage Promotion\na) Survival: Binomial(N, s_month(z, e)), β_0 downscaled\nb) Stage Promotion: Binomial(N_surv, q_month(e)), ψ_0",
           bold_first=True, palette="process")

    D.node("st3_alive", "term", 820, 495, 200, 65,
           "Active Population\nSurviving individuals\n(S, J, A per voxel)",
           bold_first=True, palette="term")

    # Annual Processes below monthly box
    # Stage 1: Reproduction
    D.node("st1_box", "lane", 335, 735, 705, 155,
           "STAGE 1: REPRODUCTION & SEED PRODUCTION", palette="loopo")
    D.node("st1_flow", "process", 360, 775, 310, 55,
           "Adult Flowering Eligibility\nBernoulli(p_flower(z_A)) onset gate")
    D.node("st1_seeds", "process", 700, 775, 315, 55,
           "Seed Production\nPoisson(N · F(z)) or floor + residual")
    D.node("st1_cost", "process", 360, 840, 655, 40,
           "Stage 4 Cost: 50% vegetative growth reduction on fruiting adults")

    # Stage 1b: 3D Wind Dispersal
    D.node("st1_disp_box", "lane", 335, 910, 705, 115,
           "STAGE 1b: 3D WIND DISPERSAL", palette="loopo")
    D.node("st1_wind", "process", 360, 950, 310, 60,
           "3D Wind Dispersal Projection\nVon Mises angle + transport (Ut, λ)")
    D.node("st1_trunc", "process", 700, 950, 315, 60,
           "Coordinate Truncation & Padding\nArray padding check (d_max, d_z)")

    # Stage 2: Seedling Establishment & Space Saturation
    D.node("st2_box", "lane", 335, 1045, 705, 160,
           "STAGE 2: SEEDLING ESTABLISHMENT & SPACE SATURATION", palette="obs")
    D.node("st2_niche", "process", 360, 1085, 310, 55,
           "Multi-Axis Climate Niche Match\nm_sp in [0, 1] suitability match")
    D.node("st2_myco", "process", 700, 1085, 315, 55,
           "Mycorrhizal Germination Filter\nBinomial(seeds, p_germ) gate")
    D.node("st2_space", "process", 360, 1150, 310, 45,
           "Space Check: C_free >= w_S")
    D.node("st2_clamp", "process", 700, 1150, 315, 45,
           "Clamping to free voxel K_xyz")

    # Stage 4: Annual Synchronization & Blending
    D.node("st4_sync", "state", 335, 1225, 705, 50,
           "STAGE 4: ANNUAL BOUNDARY SYNCHRONIZATION\nSeedling cohort size blending z_S(t+1) & update carry-over population matrices",
           bold_first=True, palette="state")

    # 4. Outputs & Evaluation
    D.node("p_out", "lane", 50, 1425, 1020, 155, "4. OUTPUTS & VIRTUAL ECOLOGIST", bold_first=True, palette="term")
    D.node("out_abund", "term", 80, 1475, 290, 80,
           "Simulated Abundance (Output)\n50-year population trajectories\n& 3D vertical profiles",
           bold_first=True, palette="term")
    D.node("out_ve", "process", 410, 1475, 295, 80,
           "Virtual Ecologist Module\nGround vs canopy sampling protocols\nHeight-dependent detection bias",
           bold_first=True, palette="process")
    D.node("out_obs", "process", 740, 1475, 300, 80,
           "Simulated Observation & Evaluation\nValidation against field records\nKS vertical displacement tests",
           bold_first=True, palette="grey")

    # CONNECTORS
    # Initialization
    D.edge("inp_canopy", "init_trees")
    D.edge("init_burnin", "st3_growth", via=[(807, 360), (575, 360)])

    # Monthly Demographics loop
    D.edge("st3_growth", "st3_surv")
    D.edge("st3_surv", "st3_alive", src_side="right", dst_side="bottom", via=[(800, 570), (920, 570)])
    D.edge("st3_alive", "st3_growth", src_side="top", dst_side="top", via=[(920, 440), (575, 440)],
           label="Month m + 1 (m < 12)", label_pos=0.5, color="#E67E22", lw=1.2)

    # From Alive adults to Reproduction
    D.edge("st3_alive", "st1_flow", src_side="bottom", dst_side="top", via=[(920, 715), (515, 715)],
           label="Month m = 12: Adult seed rain", label_pos=0.5, lw=1.2)

    # Reproduction -> Dispersal -> Establishment
    D.edge("st1_flow", "st1_seeds", orthogonal=False)
    D.edge("st1_seeds", "st1_wind", src_side="bottom", dst_side="top", via=[(857, 895), (515, 895)])
    D.edge("st1_wind", "st1_trunc", orthogonal=False)
    D.edge("st1_trunc", "st2_niche", src_side="bottom", dst_side="top", via=[(857, 1030), (515, 1030)])
    D.edge("st2_niche", "st2_myco", orthogonal=False)
    D.edge("st2_space", "st2_clamp", orthogonal=False)
    D.edge("st2_clamp", "st4_sync", src_side="bottom", dst_side="top", via=[(857, 1210), (687, 1210)])

    # Blue feedback loop: new individuals -> next generation
    D.edge("st4_sync", "st3_growth", src_side="left", dst_side="left",
           via=[(315, 1250), (315, 487)],
           label="new individuals → next generation (Year t + 1)",
           label_pos=0.45, color="#1976D2", lw=1.5)

    # Environmental input connectors
    D.edge("env_hab", "env_micro", color="#2E7D32")
    D.edge("env_micro", "st3_surv", src_side="right", dst_side="left", color="#2E7D32", lw=1.1)
    D.edge("env_nutr", "st1_flow", src_side="right", dst_side="left", via=[(310, 781), (310, 802)], color="#2E7D32", lw=1.1)
    D.edge("env_conn", "st1_wind", src_side="right", dst_side="left", color="#2E7D32", lw=1.1)
    D.edge("env_micro", "st2_niche", src_side="bottom", dst_side="left", via=[(190, 1112)], color="#2E7D32", lw=1.1)

    # Outputs & Virtual Ecologist
    D.edge("st4_sync", "out_abund", src_side="bottom", dst_side="top", via=[(687, 1315), (225, 1315)],
           label="Year t = 50 reached", label_pos=0.5, color="#B85464", lw=1.2)
    D.edge("out_abund", "out_ve", lw=1.2)
    D.edge("out_ve", "out_obs", lw=1.2)

    return D


def lookup():
    D = Diagram(760, 700, font_px=13.5)
    D.node("a", "term", 290, 15, 180, 46, "height, microenv", bold_first=False)
    D.node("b", "process", 255, 90, 250, 46, "load_height(microenv, height)", bold_first=False)
    D.node("c", "decision", 290, 165, 180, 96, "height data found?", bold_first=False)
    D.node("d", "term", 560, 190, 120, 46, "NULL", bold_first=False)
    D.node("e", "process", 220, 290, 320, 62, "Build hourly data frame\ntemperature, relative humidity, wind speed; shortwave = direct + diffuse")
    D.node("f", "decision", 270, 380, 220, 110, "microenv$.weather available?", bold_first=False)
    D.node("g1", "process", 30, 520, 300, 66, "Precipitation by calendar month\nmean per month of obs_time; wind direction: annual mean")
    D.node("g2", "process", 430, 520, 300, 66, "Precipitation: one site-wide value\nsite-wide mean, or NA if no weather record exists")
    D.node("h", "term", 270, 625, 220, 52, "Per-hour climate data frame", bold_first=False)
    D.edge("a", "b"); D.edge("b", "c")
    D.edge("c", "d", label="no", src_side="right", dst_side="left")
    D.edge("c", "e", label="yes", src_side="bottom", dst_side="top")
    D.edge("e", "f")
    D.edge("f", "g1", label="yes", src_side="left", dst_side="top", via=[(180, 435)])
    D.edge("f", "g2", label="no", src_side="right", dst_side="top", via=[(580, 435)])
    D.edge("g1", "h", src_side="bottom", dst_side="left", via=[(180, 651)])
    D.edge("g2", "h", src_side="bottom", dst_side="right", via=[(580, 651)])
    return D


def main():
    os.makedirs(OUT, exist_ok=True)
    # Output to both report/figures and manuscript for paper/figures
    targets = [
        OUT,
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "manuscript for paper", "figures")
    ]

    print("Building Program Execution Flow Diagram (ODD Scheduling Sequence)...")
    d_odd = odd_execution_flow()
    for tgt in targets:
        os.makedirs(tgt, exist_ok=True)
        d_odd.write_drawio(os.path.join(tgt, "fig2_seres_sketch.drawio"), name="ODD-Execution-Flow")
        d_odd.render_png(os.path.join(tgt, "fig2_seres_sketch.png"), dpi=260)
    print("  -> fig2_seres_sketch generated successfully!")

    print("Building Biological Causal Diagram...")
    d_causal = biological_causal_diagram()
    for tgt in targets:
        d_causal.write_drawio(os.path.join(tgt, "fig_causal_diagram.drawio"), name="Biological-Causal-Diagram")
        d_causal.render_png(os.path.join(tgt, "fig_causal_diagram.png"), dpi=260)
    print("  -> fig_causal_diagram generated successfully!")

    print("Building Climate Lookup Diagram...")
    d_lookup = lookup()
    for tgt in targets:
        d_lookup.write_drawio(os.path.join(tgt, "fig_lookup_climate.drawio"), name="Lookup-Climate")
        d_lookup.render_png(os.path.join(tgt, "fig_lookup_climate.png"), dpi=260)
    print("  -> fig_lookup_climate generated successfully!")

if __name__ == "__main__":
    main()
