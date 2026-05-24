"""
Canonical marker gene lists for Slide-tags / MERFISH mouse brain data.

MARKER_CATEGORIES drives both gene ordering and panel assignment in the dotplot.
Union of curio + merfish lists; hippocampal subregion markers appended to Glutamatergic.
"""

MARKER_CATEGORIES = {
    'Glutamatergic': [
        # ── Pan-glutamatergic ──────────────────────────────────────────────
        'Slc17a7',
        # ── Cortical L2/3 ─────────────────────────────────────────────────
        'Cux2',     # L2/3 IT
        # ── Cortical L5 ───────────────────────────────────────────────────
        'Tshz2',    # L5
        # ── Cortical L6 ───────────────────────────────────────────────────
        'Syt6',     # L6 / MH
        'Tle4',     # L6 CT
        # ── Medial habenula ───────────────────────────────────────────────
        'Tac2',     # MH Tac2 Glut
        # ── Hippocampus ───────────────────────────────────────────────────
        'Wfs1',     # CA1
        'Amigo2',   # CA2
        'Cacng5',   # CA2
        'Kcnb2',    # CA3
        'Prox1',    # DG granule cells
        'Calb1',    # DG granule / CA2
    ],
    'GABAergic': [
        # ── Pan-GABA ──────────────────────────────────────────────────────
        'Gad2',
        # ── CGE-derived ───────────────────────────────────────────────────
        'Vip',
        'Lamp5',
        'Reln',     # CGE-related
        'Calb2',    # CGE Vip/Lamp5 clusters
        # ── MGE-derived ───────────────────────────────────────────────────
        'Lhx6',     # MGE lineage factor
        'Pvalb',
        'Sst',
        # ── Striatal MSNs ─────────────────────────────────────────────────
        'Drd1',     # D1 MSNs
        'Drd2',     # D2 MSNs
    ],
    'Non-neuronal': [
        # ── Astrocytes ────────────────────────────────────────────────────
        'Slc1a3',
        # ── Ependymal ─────────────────────────────────────────────────────
        'Foxj1',
        # ── OPC ───────────────────────────────────────────────────────────
        'Pdgfra',
        'Cspg4',
        # ── Oligodendrocytes ──────────────────────────────────────────────
        'Plp1',
        # ── VLMC ──────────────────────────────────────────────────────────
        'Igf2',
        # ── Pericytes ─────────────────────────────────────────────────────
        'Vtn',      # Curio
        'Pdgfrb',   # MERFISH
        'Anpep',    # MERFISH
        # ── Endothelial ───────────────────────────────────────────────────
        'Flt1',
        'Cldn5',    # Curio
        'Pecam1',   # MERFISH
        # ── Microglia ─────────────────────────────────────────────────────
        'Cx3cr1',
    ],
}


def get_marker_union():
    """Deduplicated union across all categories, preserving order."""
    seen = set()
    result = []
    for genes in MARKER_CATEGORIES.values():
        for g in genes:
            if g not in seen:
                seen.add(g)
                result.append(g)
    return result


def get_gene_panel_map():
    """Returns {gene: broad_class} for all genes in MARKER_CATEGORIES."""
    mapping = {}
    for bc, genes in MARKER_CATEGORIES.items():
        for g in genes:
            if g not in mapping:
                mapping[g] = bc
    return mapping
