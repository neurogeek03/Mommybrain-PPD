"""
UCell signature scoring for CORT-altered glial states.

Literature-curated signatures (one per biological axis) scored on all glial
subtypes. Signature definitions and per-cell-type applicability live in
docs/gene_signatures.md — interpret each score in the context of the
expected target cell type (e.g. Myelin_identity is meaningful in oligos;
OPC_proliferation in OPCs; OXPHOS in endothelium).

Output: data/ucell_scores_v3.csv  (long format: barcode, cell_type, treatment, scores)

Run with:
    uv run python code/01_score_ucell.py
"""

import os
os.environ["NUMBA_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import anndata as ad
import pandas as pd
from pyucell import compute_ucell_scores

H5AD = "data/slide_tags_129493.h5ad"
OUT_CSV = "data/ucell_scores_v4.csv"

CELL_TYPE_COL = "subclass_name"
TREATMENT_COL = "treatment"

# 5 non-neuronal + 2 neuronal populations.
# DG granule cells (hippocampal) and Pvalb+ interneurons added to test
# the same axes in neurons. Note: cell-type-specific signatures
# (Myelin_identity, OPC_*, Astro_*, BBB_*) will be near-zero in neurons.
CELL_TYPES = [
    "327 Oligo NN",
    "326 OPC NN",
    "318 Astro-NT NN",
    "319 Astro-TE NN",
    "333 Endo NN",
    "037 DG Glut",     # dentate gyrus granule cells
    "052 Pvalb Gaba",  # parvalbumin+ interneurons
]

# ─────────────────────────────────────────────────────────────────────────────
# Signatures — source: docs/gene_signatures.md
# Comments mark expected direction in CORT and primary target cell type(s).
# DOWN signatures: low score under CORT == loss of that program.
# ─────────────────────────────────────────────────────────────────────────────
SIGNATURES = {
    # ── A. Stress / Glucocorticoid axes ──────────────────────────────────────
    # § A1 — direct GRE-containing GR targets. UP in CORT in all 5 cell types.
    # Note: rat dataset — Dusp1 not annotated, using Dusp5/6 (also GR-inducible MAPK phosphatases).
    "GR_activation": [
        "Fkbp5", "Tsc22d3", "Sgk1", "Per1", "Klf9", "Ddit4", "Dusp5", "Dusp6",  # Tier 1
        "Hif3a", "Zbtb16", "Errfi1", "Bcl6", "Nfkbia",                            # Tier 2
    ],
    # § A2 — anti-inflammatory GR brake on NF-κB / JAK-STAT. UP, strongest in Endo + Astro-TE.
    "GR_brake": [
        "Tsc22d3", "Nfkbia", "Bcl6", "Socs1", "Dusp5", "Dusp6",
        "Foxo3", "Zfp36", "Klf2",
    ],
    # § A3 — integrated stress response + hypoxia coupling. UP in all 5.
    "ISR_hypoxia": [
        "Atf4", "Atf5", "Ddit3", "Trib3",
        "Ddit4", "Deptor",
        "Eif2ak3", "Eif2ak4",
        "Hspa5", "Ppp1r15a",
        "Hif3a", "Errfi1",
    ],

    # ── B. Metabolic axes ────────────────────────────────────────────────────
    # § B1 — glucose-off → FAO + BCAA catabolism ON. UP in Oligo primarily.
    "Metabolic_reprogramming": [
        # glucose oxidation off
        "Pdk4", "Pdk2",
        # lipolysis
        "Pnpla2", "Lipe", "Plin2",
        # mitochondrial FAO
        "Cpt1a", "Cpt2", "Acadm", "Acadl", "Acadvl", "Hadha", "Hadhb",
        # peroxisomal FAO (myelin-lipid turnover)
        "Acox1", "Hsd17b4", "Ehhadh", "Slc25a17",
        # FAO regulators
        "Ppara", "Ppargc1a", "Mlycd",
        # BCAA catabolism (brain-restricted Bcat1)
        "Bcat1", "Bcat2", "Bckdha", "Bckdhb", "Dbt", "Gpt2", "Got1", "Got2",
    ],
    # § B1 (DOWN) — anabolic lipid / myelin-lipid synthesis. DOWN in CORT-Oligo.
    # Note: Scd1 not in this rat annotation; Scd2 retained.
    "Lipid_anabolism_loss": [
        "Srebf1", "Srebf2", "Fasn", "Acaca", "Acacb",
        "Elovl1", "Elovl6", "Scd2",
        "Hmgcr", "Hmgcs1", "Sqle",
        "Aspa",
    ],
    # § B2 — mitochondrial OXPHOS Complex I/II/III/IV/V + carriers. DOWN in Endo.
    # Note: Complex V uses updated HGNC symbols (Atp5f1a/b/c/d/e, Atp5mc1, Atp5po, Atp5pb).
    # Cox8a not annotated in this rat dataset; dropped.
    "OXPHOS": [
        # Complex I
        "Ndufa1", "Ndufa2", "Ndufa9", "Ndufa10", "Ndufa12",
        "Ndufb2", "Ndufb5", "Ndufb8", "Ndufb10",
        "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs7", "Ndufs8",
        "Ndufv1", "Ndufv2",
        # Complex II
        "Sdha", "Sdhb", "Sdhc", "Sdhd",
        # Complex III
        "Uqcrb", "Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrq", "Cyc1",
        # Complex IV
        "Cox4i1", "Cox5a", "Cox5b", "Cox6a1", "Cox6b1", "Cox6c", "Cox7a2", "Cox7b",
        # Complex V (new HGNC nomenclature)
        "Atp5f1a", "Atp5f1b", "Atp5f1c", "Atp5f1d", "Atp5f1e",
        "Atp5mc1", "Atp5po", "Atp5pb", "Atp5pd", "Atp5pf",
        # Mitochondrial carriers
        "Slc25a4", "Slc25a5", "Vdac1", "Vdac2", "Vdac3",
    ],

    # ── C. Cell-type identity axes (use as DOWN — loss under CORT) ───────────
    # § C1 — mature myelinating oligodendrocyte identity. DOWN in CORT-Oligo.
    # Note: rat annotation uses Ugt8 (no -a suffix).
    "Myelin_identity": [
        "Mbp", "Plp1", "Mog", "Mag", "Mobp", "Mal",
        "Opalin",
        "Cnp", "Ugt8",
        "Elovl6", "Elovl7", "Fa2h",
        "Aspa",
        "Sox10", "Myrf", "Olig2", "Nkx2-2",
    ],
    # § C2 — OPC proliferation / cell cycle. DOWN in CORT-OPC.
    "OPC_proliferation": [
        # S-phase / replication
        "Mki67", "Top2a",
        "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6", "Mcm7",
        "Pcna", "Dut", "Tyms", "Rrm1", "Rrm2", "Fancd2",
        # G2/M
        "Tpx2", "Cenpf", "Cenpa", "Cenpe",
        "Ccnb1", "Ccnb2", "Cdk1",
        "Birc5", "Aurkb", "Diaph3", "Atad2", "Dhps",
    ],
    # § C3 — pan-reactive astrocyte (Liddelow 2017 / Zamanian 2012). UP in CORT astros.
    "Astro_pan_reactive": [
        "Gfap", "Vim", "Lcn2",
        "Steap4", "S1pr3", "Cd44",
        "Osmr", "Cp", "Serpina3n",
        "Aspg", "Timp1", "Tspo", "Hsbp1",
    ],
    # § C3 — A1 neurotoxic astrocyte. UP if CORT drives an A1 skew.
    # Note: mouse-specific H2-T23/H2-D1 (MHC-I) have many rat RT1-* orthologs (1:many mapping
    # not safe to substitute); Iigp1 and Ugt1a1 also lack 1:1 rat orthologs in this annotation.
    # Remaining markers are core canonical A1.
    "Astro_A1_neurotoxic": [
        "C3", "Serping1", "Gbp2",
        "Ggta1", "Srgn", "Fbln5",
        "Amigo2", "Psmb8",
    ],
    # § C3 — A2 neuroprotective astrocyte. UP if CORT drives an A2 skew.
    "Astro_A2_protective": [
        "Clcf1", "Tgm1", "Ptx3",
        "S100a10", "Sphk1", "Cd109",
        "Ptgs2", "Emp1",
        "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14",
    ],
    # § C4 — BBB efflux + transport machinery. UP in CORT-Endo.
    "BBB_transport": [
        "Abcg2", "Abcb1a", "Abcb1b", "Abcc4",      # efflux
        "Slc2a1", "Slc7a5", "Slc16a1", "Slc19a3",  # solute carriers
        "Mfsd2a", "Lpl", "Timp3",
    ],
    # § C4 — vascular identity TFs / receptors. DOWN under prolonged CORT.
    "BBB_vascular_identity": [
        "Kdr", "Lama4", "Pdgfb", "Thra", "Notch1", "Notch4",
    ],

    # ── D. Neurotransmission / synapse axes (astrocyte-specific, DOWN) ───────
    # § D1 — astrocyte glutamatergic synapse support. DOWN in CORT-Astro-TE.
    "Astro_glutamate_support": [
        "Slc1a2", "Slc1a3", "Glul",
        "Gpc4", "Gpc6", "Sparcl1", "Sparc",
        "Thbs1", "Thbs2", "Chrdl1",
        "Grm3", "Grm5",
        "Frrs1l",
        "Epha4", "Ephb1",
    ],
    # § D2 — astrocyte GABAergic support. DOWN in CORT-Astro-NT (PPD-relevant).
    "Astro_GABA_support": [
        "Gabbr1", "Gabbr2",
        "Gabra2", "Gabra4", "Gabrb1", "Gabrg2",
        "Slc6a1", "Slc6a11", "Slc6a13",
        "Abat", "Aldh5a1",
        "Adcy8", "Best1",
    ],

    # ── E. Cross-cutting modules ─────────────────────────────────────────────
    # § E1 — shared glial cytoskeletal remodeling. UP in OPC, both Astros.
    "Glial_cyto_remodel": [
        "Tmod1", "Tmod2", "Mical2", "Svil", "Filip1",
        "Pak1", "Cyria", "Asap2", "Minar1",
    ],
    # § E2 — core circadian clock. DOWN in CORT-Endo (Per1 caveat: direct GR target).
    # Note: rat annotation uses Bmal1 (not Arntl).
    "Clock_genes": [
        "Bmal1", "Clock", "Npas2",
        "Per1", "Per2", "Per3",
        "Cry1", "Cry2",
        "Nr1d1", "Nr1d2",
        "Rora", "Rorb", "Rorc",
        "Dbp", "Tef", "Hlf",
    ],
    # § E3 — adenosine activity-sensing / suppression. UP in parenchymal glia.
    "Adenosine_activity_sensing": [
        "Adora1", "Adora2a", "Adora2b", "Adora3",
        "Nt5e", "Ada", "Entpd1",
    ],
}


def main():
    print("Loading AnnData...")
    adata = ad.read_h5ad(H5AD)
    print(f"  Full dataset: {adata.shape}")

    adata.var_names = adata.var["gene_symbol"].astype(str)
    adata.var_names_make_unique()
    adata.var.index.name = None

    # Verify signature genes are present
    for sig_name, genes in SIGNATURES.items():
        missing = [g for g in genes if g not in adata.var_names]
        found = [g for g in genes if g in adata.var_names]
        print(f"\n{sig_name}: {len(found)}/{len(genes)} genes found"
              + (f"  MISSING: {missing}" if missing else ""))

    # Subset to selected cell types and score
    mask = adata.obs[CELL_TYPE_COL].isin(CELL_TYPES)
    sub = adata[mask].copy()
    print(f"\nScoring {sub.n_obs} cells across {len(CELL_TYPES)} subtypes "
          f"x {len(SIGNATURES)} signatures...")

    compute_ucell_scores(sub, SIGNATURES, n_jobs=1)

    score_cols = [f"{s}_UCell" for s in SIGNATURES]

    # Build output: barcode, cell_type, treatment, scores
    out = sub.obs[[CELL_TYPE_COL, TREATMENT_COL]].copy()
    out.index.name = "barcode"
    for col in score_cols:
        out[col] = sub.obs[col].values

    out.to_csv(OUT_CSV)
    print(f"\nSaved {OUT_CSV}  ({len(out)} rows)")
    print(out.groupby([CELL_TYPE_COL, TREATMENT_COL])[score_cols].mean().round(4).to_string())


if __name__ == "__main__":
    main()
