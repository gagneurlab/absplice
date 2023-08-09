#!/usr/bin/env python3
# requires "interpret==0.2.7" "interpret_core==0.2.7" "ebm2onnx==1.3"
import numpy as np
import onnx
import onnx.helper
import ebm2onnx

import pickle
import json

with open("absplice/precomputed/AbSplice_DNA.pkl", "rb") as fd:
    absplice_dna_model = pickle.load(fd)

print(json.dumps(dict(zip(absplice_dna_model.feature_names, absplice_dna_model.feature_types)), indent=2))

onnx_model = ebm2onnx.to_onnx(
    absplice_dna_model,
    {
        'delta_logit_psi': 'double',
        'delta_psi': 'double',
        'delta_score': 'double',
        'splice_site_is_expressed': 'int',
        # 'delta_psi_cat': 'double'
    },
    predict_proba=True
)

onnx.save_model(onnx_model, 'absplice/precomputed/AbSplice_DNA.onnx')

with open("absplice/precomputed/AbSplice_RNA.pkl", "rb") as fd:
    absplice_rna_model = pickle.load(fd)

print(json.dumps(dict(zip(absplice_rna_model.feature_names, absplice_rna_model.feature_types)), indent=2))

onnx_model = ebm2onnx.to_onnx(
    absplice_rna_model,
    {
        'delta_logit_psi': 'double',
        'delta_psi': 'double',
        'delta_psi_cat': 'double',
        'delta_score': 'double',
        'pValueGene_g_minus_log10': "double",
        'splice_site_is_expressed': 'int',
    },
    predict_proba=True
)

onnx.save_model(onnx_model, 'absplice/precomputed/AbSplice_RNA.onnx')

with open("absplice/precomputed/ABSPLICE_DNA_with_CADD_Splice.pkl", "rb") as fd:
    absplice_dna_cadd_model = pickle.load(fd)

print(json.dumps(dict(zip(absplice_dna_cadd_model.feature_names, absplice_dna_cadd_model.feature_types)), indent=2))

onnx_model = ebm2onnx.to_onnx(
    absplice_dna_cadd_model,
    {
        'PHRED': 'double',
        'delta_logit_psi': 'double',
        'delta_psi': 'double',
        'delta_score': 'double',
        'splice_site_is_expressed': 'int',
    },
    predict_proba=True
)
onnx.save_model(onnx_model, 'absplice/precomputed/ABSPLICE_DNA_with_CADD_Splice.onnx')
