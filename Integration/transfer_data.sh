#!/bin/bash
# Run this from Trillium to copy data files to Nibi
# Usage: bash transfer_data.sh

DEST="mfafouti@nibi.alliancecan.ca:/scratch/mfafouti/Mommybrain/Integration/data/"

## h5ad files
rsync -avP \
  "/scratch/mfafouti/Mommybrain/Slide_tags/Filtering/NEW_list_merged_filtered/Merged_n=207903_all_metadata_slide_tags.h5ad" \
  "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/out/semi_filtered_neurons/adata_neuron_subset_181663_singlets_score_330.h5ad" \
  "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/out/adata_nn_subset_78689_singlets_score_330.h5ad" \
  "$DEST"

## Slide-tags spatial coordinates
rsync -avP \
  /scratch/mfafouti/Mommybrain/Slide_tags/Spatial/coords_BC3.csv \
  /scratch/mfafouti/Mommybrain/Slide_tags/Spatial/coords_BC9.csv \
  /scratch/mfafouti/Mommybrain/Slide_tags/Spatial/coords_BC13.csv \
  /scratch/mfafouti/Mommybrain/Slide_tags/Spatial/coords_BC14.csv \
  /scratch/mfafouti/Mommybrain/Slide_tags/Spatial/coords_BC15.csv \
  /scratch/mfafouti/Mommybrain/Slide_tags/Spatial/coords_BC28.csv \
  "$DEST"
