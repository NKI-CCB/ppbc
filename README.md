# Post-partum breast cancer

This repository contains the code for the analysis of RNA-seq data from the KWF project "Postpartum breast cancer diagnosed during involution: a distinct entity with unique clinicopathological, molecular and immunological features.

## Project intro/objectives

From the KWF project proposal:

"Postpartum breast cancer (PPBC) represents an under-recognised yet high-risk subset of young women with breast cancer. Breast cancers that occur 1 to 2 years postpartum are associated with worse survival rates and a 2-fold increased risk for metastasis compared with breast cancers diagnosed in young, premenopausal women during or outside pregnancy.

The poorer prognosis and increased risk for metastasis and death of PPBC cannot solely be explained by altered clinical risk factors, like delayed diagnosis and deregulation of gestational hormones during or shortly after pregnancy. Although the exact molecular mechanisms remain unclear, it has been hypothesized that changes in the breast, uniquely occurring in the postpartum window, are underlying the poor prognosis of postpartum breast cancer. During pregnancy, the mammary gland epithelium undergoes remarkable proliferation and differentiation in order to prepare for lactation. After delivery, or (if breastfeeding is given) after lactation, the fully differentiated gland regresses to its pre-pregnant state, a process called involution. In rodents, the process of mammary gland involution shares striking similarities with the tissue-remodelling programs that are activated during wound healing and inflammation. Rodent mammary gland involution results in numerous changes in extracellular matrix (ECM) organization, which in turn stimulate the influx and activation of immunosuppressive leukocytes, like regulatory T-cells and tumour-associated macrophages with T-cell-suppressive function, and promote tumour cell escape from the mammary gland. 

In women, the correlation of immune cell infiltration and poor prognosis of breast cancer in general has been well studied, and has been driving investigations into immunotherapy as a means of breast cancer treatment. However, extensive characterization of the (immuno)biology of the mammary gland of postpartum breast cancer
patients, is lacking. 

In this research project we AIM to investigate whether PPBC, specifically diagnosed during mammary gland involution (PPBC-inv), has worse survival and metastatic rates and exhibits unique molecular and immunological characteristics when compared to breast cancers diagnosed during pregnancy (PrBC) or outside the context of a pregnancy (non-PrBC, nulliparous women). In a supplementary exploratory analysis we will examine how features of PPBC-lac (i.e. postpartum breast cancers diagnosed during lactation) relate to those of PrBC and PPBC-inv, the first two types of breast cancer being diagnosed in pre-involutional women."

Preliminary findings have shown that cancers diagnosed during involution do indeed have substantially worse outcomes than those diagnosed during pregnancy or in nulliparous women.

## Research questions

* What is the molecular mechanism behind the poorer outcome in involuting women?
* Is there an immune signature that is predictive of clinical outcome?
* Does (extended) breastfeeding confer better outcome in PPBC?
* Does the duration of involution (between cessation of breastfeeding and diagnosis) impact clinical outcome?

## Approach

This project utilizes a three-pronged approach:

1) Survival analysis based on clinical outcomes within the PPBC cohort
2) RNAseq-based analyses on FFPE-preserved primary tumors from a subset of patients within the PPBC cohort
3) Spatial analysis from Vectra multiplex panels on a subset of the aforementioned patients

## Directory structure

Based on [cookiecutter data science](https://drivendata.github.io/cookiecutter-data-science/) but does not perfectly conform to that formula.

WIP (fill in dir structure when final)

RNAseq reports are in the parent `reports` directory. Spatial reports are in a sub-folder: `reports/spatial`.
