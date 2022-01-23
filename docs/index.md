# CancerVar & OPAI

## CancerVar

CancerVar is a bioinformatics software tool for clinical interpretation of somatic variants by the AMP/ASCO/CAP 2017 guidelines. The input to CancerVar is an annotated file generated from ANNOVAR, while the output of CancerVar is the classification of variants into "Tier_I_strong","Tier_II_potential","Tier_IV_benign" and "Tier_III_Uncertain", together with detailed evidence code. CancerVar can have automated interpretation and manual interpretation.

- **Automated interpretation**: The standards and guidelines for the clinical interpretation of sequence variants in cancer, based on 12 criteria. The CancerVar can generate automated interpretation on 10 criteria: (CBP1,CBP2,CBP3,CBP4,CBP7,CBP8,CBP9,CBP10,CBP11,CBP12)

- **Manual interpretation**: The rest of  criteria (CBP5, CBP6) requires user input in the manual adjustment step.  The user can provide these criteria or other criteria in the evidence file.

## OPAI

OPAI is a python script for oncogenic prioritization by artificial intelligence (OPAI). Using deep learning-based approach, OPAI provides probability score  to determine the oncogenicity of a variant, using 12 clinical evidence features from CancerVar and 19 functional features with scoring metrics predicted by computational tools.

## Pipeline

 Check [here](misc/whatsnew.md) to see what is new in CancerVar & OPAI.

---

![new](img/new.png) 2018Nov: The web server of CancerVar is  online now [here](http://cancervar.wglab.org).

![new](img/new.png) 2018Jan: The GitHub repository for CancerVar is created.

---

## Reference

- Quan Li, Zilin Ren, Kajia Cao, Marilyn M. Li, Yunyun Zhou and Kai Wang. CancerVar: an Artificial Intelligence empowered platform for clinical interpretation of somatic mutations in cancer.

- Li MM, Datto M, Duncavage EJ, Kulkarni S, Lindeman NI, Roy S, Tsimberidou AM, Vnencak-Jones CL, Wolff DJ, Younes A, Nikiforova MN. Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists.

- Li Q, Wang K. InterVar: Clinical Interpretation of Genetic Variants by the 2015 ACMG-AMP Guidelines. _American Journal of Human Genetics_, 2:267-280, 2017 


