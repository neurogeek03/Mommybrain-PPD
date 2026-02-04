This file contains notes to help the Gemini agent provide focused responses for the LIANA project.

**Project Goal:** The user wants to perform cell-cell communication (CCC) analysis using the `liana` Python package.

**Key Components:**
*   **Core Tool:** `liana` (Python version).
*   **Main Script:** `scripts/CCC_liana.py` is the primary entry point for the analysis.
*   **Data Format:** The project uses AnnData (`.h5ad`) files for single-cell data, located in the `data/` directory.
*   **Dependencies:** The environment is managed by Conda via `liana_env.yml`.

**User's Current Problem (2026-02-02):**
*   A gene filtering step in `scripts/CCC_liana.py` is removing all genes for certain cell types (e.g., "318 Astro-NT NN").
*   The user observes non-zero counts for genes within these cell types, indicating the filtering logic is too strict or is not behaving as expected.
*   **Next Step:** Investigate the filtering logic in `scripts/CCC_liana.py` and potentially `scripts/functions.py` to identify the responsible code and parameters.