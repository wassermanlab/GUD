
# Genomic Universal Database

> This is a database which centralizes and unifies genomic data in a universal manner for specific reference datasets.


## Manifest
+ GUD/ORM - Object relational mapping classes
+ GUD/scripts - Scripts for upserting data into GUD
+ 

## Requirements
```
pip install Biopython
pip install mysqlclient
pip install interval-binning
pip install sqlalchemy
```

## How to append GUD as submodule
```
MAIN_PROJECT_FOLDER/
    LICENCE
    README.md
    .gitignore
    PROJECT_FOLDER/
        __init__.py
        lib/
        src/
        GUD/
```
when inside `PROJECT_FOLDER` type:
```
git submodule add --force git@github.com:oriolfornes/GUD.git
```

## Gene representation
![](https://github.com/oriolfornes/GUD/blob/master/GUD-Gene.png)

## Important Notes
+ Everything in GUD is zero-based open-bookend to keep a standard with BED format
+ Data currently within GUD found [here](https://docs.google.com/document/d/1Jjug_gvsTZk-E1L2sMa8AflUq3SW9RVPOIuNtsQjQZM/edit)



```
mysql> show tables;
+----------------------+
| Tables_in_hg19       |
+----------------------+
| chrom_size           |
| conservation         |
| dna_accessibility    |
| enhancer             |
| gene                 |
| histone_modification |
| rmsk                 |
| tad                  |
| tf_binding           |
+----------------------+
9 rows in set (0.03 sec)

+----------+

mysql> select count(*) from dna_accessibility;
+----------+
| count(*) |
+----------+
| 58432698 |
+----------+
1 row in set (0.00 sec)

mysql> select distinct(cell_or_tissue) from dna_accessibility;
+-------------------------------------------------------------------------+
| cell_or_tissue                                                          |
+-------------------------------------------------------------------------+
| neuron, hippocampus                                                     |
| glia, orbitofrontal cortex                                              |
| glia, dorsolateral prefrontal cortex                                    |
| neuron, orbitofrontal cortex                                            |
| glia, putamen                                                           |
| neuron, inferior temporal cortex                                        |
| glia, primary motor cortex                                              |
| glia, mediodorsal thalamus                                              |
| glia, inferior temporal cortex                                          |
| glia, hippocampus                                                       |
| glia, amygdala                                                          |
| glia, ventrolatelar prefrontal cortex                                   |
| glia, superior temporal cortex                                          |
| glia, primary visual cortex                                             |
| neuron, insula                                                          |
| glia, anterior cingulate cortex                                         |
| neuron, superior temporal cortex                                        |
| neuron, amygdala                                                        |
| neuron, primary visual cortex                                           |
| neuron, anterior cingulate cortex                                       |
| glia, insula                                                            |
| neuron, ventrolatelar prefrontal cortex                                 |
| neuron, dorsolateral prefrontal cortex                                  |
| neuron, putamen                                                         |
| neuron, mediodorsal thalamus                                            |
| neuron, nucleus accumbens                                               |
| neuron, primary motor cortex                                            |
| endothelial cell of umbilical vein primary cell                         |
| HeLa-S3 cell line                                                       |
| glia, nucleus accumbens                                                 |
| HepG2 cell line                                                         |
| GM12878 cell line                                                       |
| SJCRH30 cell line                                                       |
| hindlimb muscle tissue                                                  |
| T47D cell line                                                          |
| muscle of back tissue                                                   |
| frontal cortex tissue                                                   |
| HTR-8/SVneo cell line                                                   |
| CD14-positive monocyte primary cell                                     |
| large intestine tissue                                                  |
| left lung tissue                                                        |
| HuH-7.5 cell line                                                       |
| psoas muscle tissue                                                     |
| stomach tissue                                                          |
| iPS-NIHi7 induced pluripotent stem cell line                            |
| heart tissue                                                            |
| medulloblastoma cell line                                               |
| renal cortex interstitium tissue                                        |
| epidermal melanocyte primary cell                                       |
| naive B cell primary cell                                               |
| B cell primary cell                                                     |
| coronary artery tissue                                                  |
| kidney tissue                                                           |
| HPDE6-E6E7 cell line                                                    |
| osteoblast primary cell                                                 |
| pancreas tissue                                                         |
| tibial nerve tissue                                                     |
| muscle of arm tissue                                                    |
| right lung tissue                                                       |
| left renal cortex interstitium tissue                                   |
| HS-27A cell line                                                        |
| hematopoietic multipotent progenitor cell in vitro differentiated cells |
| fibroblast of gingiva primary cell                                      |
| astrocyte of the spinal cord primary cell                               |
| right renal pelvis tissue                                               |
| breast epithelium tissue                                                |
| M059J cell line                                                         |
| foreskin keratinocyte primary cell                                      |
| fibroblast of dermis primary cell                                       |
| fibroblast of mammary gland primary cell                                |
| HuH-7 cell line                                                         |
| CD1c-positive myeloid dendritic cell primary cell                       |
| astrocyte of the cerebellum primary cell                                |
| COLO829 cell line                                                       |
| GM19238 cell line                                                       |
| NB4 cell line                                                           |
| ascending aorta tissue                                                  |
| thymus tissue                                                           |
| MCF-7 cell line                                                         |
| OCI-LY7 cell line                                                       |
| ELR cell line                                                           |
| IMR-90 cell line                                                        |
| chorion tissue                                                          |
| GM10248 cell line                                                       |
| GM03348 cell line                                                       |
| GM10266 cell line                                                       |
| Karpas-422 cell line                                                    |
| left cardiac atrium tissue                                              |
| right renal cortex interstitium tissue                                  |
| T-helper 2 cell primary cell                                            |
| AG20443 cell line                                                       |
| astrocyte of the hippocampus primary cell                               |
| cerebellum tissue                                                       |
| GM19240 cell line                                                       |
| hepatic stellate cell primary cell                                      |
| left kidney tissue                                                      |
| mesenchymal stem cell in vitro differentiated cells                     |
| renal pelvis tissue                                                     |
| GM12891 cell line                                                       |
| iPS-NIHi11 induced pluripotent stem cell line                           |
| spinal cord tissue                                                      |
| T-cell primary cell                                                     |
| dedifferentiated amniotic fluid mesenchymal stem cell stem cell         |
| H1-hESC stem cell                                                       |
| Daoy cell line                                                          |
| foreskin melanocyte primary cell                                        |
| heart left ventricle tissue                                             |
| olfactory neurosphere cell line in vitro differentiated cells           |
| cardiac mesoderm in vitro differentiated cells                          |
| esophagus muscularis mucosa tissue                                      |
| GM12865 cell line                                                       |
| HEK293T cell line                                                       |
| 8988T cell line                                                         |
| H54 cell line                                                           |
| T-helper 1 cell primary cell                                            |
| astrocyte primary cell                                                  |
| HCT116 cell line                                                        |
| hepatocyte in vitro differentiated cells                                |
| lung microvascular endothelial cell primary cell                        |
| small intestine tissue                                                  |
| urothelium cell line cell line                                          |
| esophagus squamous epithelium tissue                                    |
| GM06990 cell line                                                       |
| naive thymus-derived CD4-positive, alpha-beta T cell primary cell       |
| A549 cell line                                                          |
| adrenal gland tissue                                                    |
| AG10803 cell line                                                       |
| amniotic epithelial cell primary cell                                   |
| arm bone tissue                                                         |
| BE2C cell line                                                          |
| body of pancreas tissue                                                 |
| brain tissue                                                            |
| Caco-2 cell line                                                        |
| cardiac fibroblast primary cell                                         |
| CD4-positive, alpha-beta T cell primary cell                            |
| CD8-positive, alpha-beta T cell primary cell                            |
| common myeloid progenitor, CD34-positive primary cell                   |
| dermis blood vessel endothelial cell primary cell                       |
| endodermal cell in vitro differentiated cells                           |
| eye tissue                                                              |
| femur tissue                                                            |
| fibroblast of lung primary cell                                         |
| fibroblast of the aortic adventitia primary cell                        |
| fibroblast of villous mesenchyme primary cell                           |
| gastrocnemius medialis tissue                                           |
| GM23248 cell line                                                       |
| GM23338 induced pluripotent stem cell line                              |
| H7-hESC stem cell                                                       |
| H9 stem cell                                                            |
| HAP-1 cell line                                                         |
| heart right ventricle tissue                                            |
| islet precursor cell in vitro differentiated cells                      |
| Jurkat clone E61 cell line                                              |
| leg bone tissue                                                         |
| liver tissue                                                            |
| LNCaP clone FGC cell line                                               |
| lower leg skin tissue                                                   |
| lung tissue                                                             |
| mesendoderm in vitro differentiated cells                               |
| muscle of leg tissue                                                    |
| muscle of trunk tissue                                                  |
| natural killer cell primary cell                                        |
| neural progenitor cell in vitro differentiated cells                    |
| omental fat pad tissue                                                  |
| ovary tissue                                                            |
| Peyer's patch tissue                                                    |
| placenta tissue                                                         |
| prostate gland tissue                                                   |
| retina tissue                                                           |
| right atrium auricular region tissue                                    |
| right kidney tissue                                                     |
| right lobe of liver tissue                                              |
| sigmoid colon tissue                                                    |
| SK-MEL-5 cell line                                                      |
| skeletal muscle myoblast primary cell                                   |
| spleen tissue                                                           |
| testis tissue                                                           |
| thyroid gland tissue                                                    |
| tibial artery tissue                                                    |
| tongue tissue                                                           |
| transverse colon tissue                                                 |
| trophoblast cell primary cell                                           |
| upper lobe of left lung tissue                                          |
| urinary bladder tissue                                                  |
| uterus tissue                                                           |
| vagina tissue                                                           |
| cardiac muscle cell primary cell                                        |
| fibroblast of peridontal ligament primary cell                          |
| occipital lobe tissue                                                   |
| iPS DF 19.11 induced pluripotent stem cell line                         |
| T-helper 17 cell primary cell                                           |
| EL cell line                                                            |
| mammary epithelial cell primary cell                                    |
| L1-S8R induced pluripotent stem cell line                               |
| trophoblast cell in vitro differentiated cells                          |
| fibroblast of pulmonary artery primary cell                             |
| CMK cell line                                                           |
| WERI-Rb-1 cell line                                                     |
| foreskin fibroblast primary cell                                        |
| regulatory T cell primary cell                                          |
| colon tissue                                                            |
| endometrium tissue                                                      |
| RCC 7860 cell line                                                      |
| epithelial cell of esophagus primary cell                               |
| GM18507 cell line                                                       |
| G401.6 cell line                                                        |
| A204.1 cell line                                                        |
| TTC549 cell line                                                        |
| skin fibroblast primary cell                                            |
| K562 cell line                                                          |
| PC-9 cell line                                                          |
| CWRU1 induced pluripotent stem cell line                                |
| AG08395 cell line                                                       |
| myotube in vitro differentiated cells                                   |
| HK-2 cell line                                                          |
| renal cortical epithelial cell primary cell                             |
| AG08396 cell line                                                       |
| bronchial epithelial cell primary cell                                  |
| choroid plexus epithelial cell primary cell                             |
| iris pigment epithelial cell primary cell                               |
| KBM-7 cell line                                                         |
| hepatocyte primary cell                                                 |
| brain pericyte primary cell                                             |
| WI38 cell line                                                          |
| HL-60 cell line                                                         |
| GM04503 cell line                                                       |
| forelimb muscle tissue                                                  |
| non-pigmented ciliary epithelial cell primary cell                      |
| iPS DF 6.9 induced pluripotent stem cell line                           |
| RWPE1 cell line                                                         |
| AG04449 cell line                                                       |
| NT2/D1 cell line                                                        |
| left renal pelvis tissue                                                |
| iPS DF 4.7 induced pluripotent stem cell line                           |
| keratinocyte primary cell                                               |
| retinal pigment epithelial cell primary cell                            |
| skeletal muscle cell primary cell                                       |
| HFF-Myc cell line                                                       |
| AG04450 cell line                                                       |
| AG09309 cell line                                                       |
| AG09319 cell line                                                       |
| Panc1 cell line                                                         |
| neural stem progenitor cell in vitro differentiated cells               |
| smooth muscle cell of the brain vasculature primary cell                |
| BJ cell line                                                            |
| brain microvascular endothelial cell primary cell                       |
| dermis microvascular lymphatic vessel endothelial cell primary cell     |
| fibroblast of the conjunctiva primary cell                              |
| germinal center tissue                                                  |
| GM04504 cell line                                                       |
| HS-5 cell line                                                          |
| iPS DF 19.7 induced pluripotent stem cell line                          |
| skin of body tissue                                                     |
| stromal cell of bone marrow primary cell                                |
| kidney epithelial cell primary cell                                     |
| SK-N-MC cell line                                                       |
| GM12864 cell line                                                       |
| fibroblast of skin of abdomen primary cell                              |
| glomerular endothelial cell primary cell                                |
| GM12892 cell line                                                       |
| GM19239 cell line                                                       |
| epithelial cell of proximal tubule primary cell                         |
| pulmonary artery endothelial cell primary cell                          |
| epithelial cell of prostate primary cell                                |
| islet of Langerhans tissue                                              |
+-------------------------------------------------------------------------+
265 rows in set (9 min 43.37 sec)

mysql> select distinct(experiment_type) from dna_accessibility;
+-----------------+
| experiment_type |
+-----------------+
| ATAC-seq        |
| FAIRE-seq       |
| DNase-seq       |
+-----------------+
3 rows in set (3 min 49.07 sec)

mysql> select distinct(source_name) from dna_accessibility;
+-------------+
| source_name |
+-------------+
| BOCA        |
| ENCODE      |
+-------------+
2 rows in set (4 min 2.73 sec)

+----------+

mysql> select count(*) from enhancer;
+----------+
| count(*) |
+----------+
|   992865 |
+----------+
1 row in set (0.00 sec)

mysql> select distinct(cell_or_tissue) from enhancer;
+--------------------------------------+
| cell_or_tissue                       |
+--------------------------------------+
| Ntera2 embryonic stem cells          |
| A549                                 |
| MCF-7                                |
| Jurkat                               |
| CD4-positive T cells                 |
| AC16                                 |
| IMR90                                |
| K562                                 |
| HeLa                                 |
| HCT-116                              |
| Ramos                                |
| LNCaP                                |
| HUVEC                                |
| B cells                              |
| HAEC                                 |
| H1 embryonic stem cells              |
| GM12878                              |
| HeLa-S3                              |
| hindbrain (rhombencephalon)          |
| heart                                |
| branchial arch                       |
| eye                                  |
| facial mesenchyme                    |
| limb                                 |
| neural tube                          |
| nose                                 |
| forebrain                            |
| somite                               |
| cranial nerve                        |
| midbrain (mesencephalon)             |
| other                                |
| ear                                  |
| dorsal root ganglion                 |
| tail                                 |
| trigeminal V (ganglion, cranial)     |
| blood vessels                        |
| genital tubercle                     |
| pancreas                             |
| liver                                |
| melanocytes                          |
| mesenchyme derived from neural crest |
+--------------------------------------+
41 rows in set (7.42 sec)

mysql> select distinct(experiment_type) from enhancer;
+-----------------+
| experiment_type |
+-----------------+
| GRO-seq         |
| PRO-seq         |
| GRO-cap         |
| 5'GRO           |
| in vivo         |
+-----------------+
5 rows in set (1 min 30.13 sec)

+----------+

mysql> select count(*) from gene;
+----------+
| count(*) |
+----------+
|    68957 |
+----------+
1 row in set (0.02 sec)

+----------+

# temp number (upserting data)

mysql> select count(*) from histone_modification;
+-----------+
| count(*)  |
+-----------+
| 202566780 |
+-----------+
1 row in set (0.06 sec)

mysql> select distinct(histone_type) from histone_modification;
+--------------+
| histone_type |
+--------------+
| H3K4me3      |
| H3K4me1      |
| H3K27ac      |
| H3K36me3     |
| H3K9me3      |
| H3K23me2     |
| H3K27me3     |
| H3F3A        |
| H4K8ac       |
| H4K20me1     |
| H4K91ac      |
| H3K9ac       |
| H3K4me2      |
| H4K5ac       |
| H3K79me2     |
| H2AFZ        |
| H3K9me2      |
| H3K79me1     |
| H3K9me1      |
| H2AK9ac      |
| H3K14ac      |
| H3K56ac      |
| H2BK5ac      |
| H3K18ac      |
| H3K23ac      |
| H3K4ac       |
| H2AK5ac      |
| H2BK12ac     |
| H2BK15ac     |
| H2BK120ac    |
| H2BK20ac     |
| H3T11ph      |
+--------------+
32 rows in set (17 min 59.88 sec)

mysql> select distinct(cell_or_tissue) from histone_modification;
+------------------------------------------------------------------------+
| cell_or_tissue                                                         |
+------------------------------------------------------------------------+
| chorionic villus tissue                                                |
| chorion tissue                                                         |
| fibroblast of breast primary cell                                      |
| luminal epithelial cell of mammary gland primary cell                  |
| T-cell primary cell                                                    |
| foreskin fibroblast primary cell                                       |
| UCSF-4 stem cell                                                       |
| H1-hESC stem cell                                                      |
| foreskin melanocyte primary cell                                       |
| germinal matrix tissue                                                 |
| mammary stem cell stem cell                                            |
| GM23248 cell line                                                      |
| amnion tissue                                                          |
| GM23338 induced pluripotent stem cell line                             |
| breast epithelium tissue                                               |
| foreskin keratinocyte primary cell                                     |
| hepatocyte in vitro differentiated cells                               |
| brain tissue                                                           |
| SK-N-SH cell line                                                      |
| heart left ventricle tissue                                            |
| K562 cell line                                                         |
| cerebellum tissue                                                      |
| Peyer's patch tissue                                                   |
| coronary artery tissue                                                 |
| NCI-H929 cell line                                                     |
| adrenal gland tissue                                                   |
| ascending aorta tissue                                                 |
| body of pancreas tissue                                                |
| gastrocnemius medialis tissue                                          |
| GM08714 cell line                                                      |
| fibroblast of villous mesenchyme primary cell                          |
| CD14-positive monocyte primary cell                                    |
| AG04450 cell line                                                      |
| SU-DHL-6 cell line                                                     |
| large intestine tissue                                                 |
| HUES6 stem cell                                                        |
| Loucy cell line                                                        |
| H9 stem cell                                                           |
| iPS DF 19.11 induced pluripotent stem cell line                        |
| common myeloid progenitor, CD34-positive primary cell                  |
| NT2/D1 cell line                                                       |
| B cell primary cell                                                    |
| HeLa-S3 cell line                                                      |
| DND-41 cell line                                                       |
| endothelial cell of umbilical vein primary cell                        |
| colonic mucosa tissue                                                  |
| mesenchymal stem cell in vitro differentiated cells                    |
| kidney epithelial cell primary cell                                    |
| IMR-90 cell line                                                       |
| lung tissue                                                            |
| BJ cell line                                                           |
| HepG2 cell line                                                        |
| bronchial epithelial cell primary cell                                 |
| esophagus muscularis mucosa tissue                                     |
| caudate nucleus tissue                                                 |
| HCT116 cell line                                                       |
| Jurkat clone E61 cell line                                             |
| A673 cell line                                                         |
| endocrine pancreas tissue                                              |
| cardiac mesoderm in vitro differentiated cells                         |
| fibroblast of dermis primary cell                                      |
| MM.1S cell line                                                        |
| CD4-positive, alpha-beta memory T cell primary cell                    |
| fibroblast of lung primary cell                                        |
| CD8-positive, alpha-beta T cell primary cell                           |
| H7-hESC stem cell                                                      |
| ectodermal cell in vitro differentiated cells                          |
| PC-9 cell line                                                         |
| VCaP cell line                                                         |
| RWPE1 cell line                                                        |
| SK-N-MC cell line                                                      |
| mammary epithelial cell primary cell                                   |
| Caco-2 cell line                                                       |
| NB4 cell line                                                          |
| esophagus tissue                                                       |
| MCF-7 cell line                                                        |
| ES-I3 stem cell                                                        |
| HUES64 stem cell                                                       |
| cardiac muscle cell in vitro differentiated cells                      |
| CD4-positive, alpha-beta T cell primary cell                           |
| DOHH2 cell line                                                        |
| GM12878 cell line                                                      |
| Karpas-422 cell line                                                   |
| OCI-LY7 cell line                                                      |
| C4-2B cell line                                                        |
| KMS-11 cell line                                                       |
| Panc1 cell line                                                        |
| gastroesophageal sphincter tissue                                      |
| OCI-LY1 cell line                                                      |
| ACC112 cell line                                                       |
| esophagus squamous epithelium tissue                                   |
| liver tissue                                                           |
| fibroblast of the aortic adventitia primary cell                       |
| LNCaP clone FGC cell line                                              |
| HFF-Myc cell line                                                      |
| GM12864 cell line                                                      |
| OCI-LY3 cell line                                                      |
| HUES48 stem cell                                                       |
| PC-3 cell line                                                         |
| A549 cell line                                                         |
| GM12865 cell line                                                      |
| cardiac fibroblast primary cell                                        |
| iPS-18a induced pluripotent stem cell line                             |
| RWPE2 cell line                                                        |
| epithelial cell of proximal tubule primary cell                        |
| cardiac muscle cell primary cell                                       |
| GM06990 cell line                                                      |
| brain microvascular endothelial cell primary cell                      |
| HL-60 cell line                                                        |
| AG09319 cell line                                                      |
| keratinocyte primary cell                                              |
| astrocyte of the spinal cord primary cell                              |
| mesendoderm in vitro differentiated cells                              |
| iPS DF 6.9 induced pluripotent stem cell line                          |
| 22Rv1 cell line                                                        |
| astrocyte of the cerebellum primary cell                               |
| AG09309 cell line                                                      |
| endodermal cell in vitro differentiated cells                          |
| epithelial cell of prostate primary cell                               |
| KOPT-K1 cell line                                                      |
| CD8-positive, alpha-beta memory T cell primary cell                    |
| layer of hippocampus tissue                                            |
| iPS-20b induced pluripotent stem cell line                             |
| kidney tissue                                                          |
| astrocyte primary cell                                                 |
| CD4-positive, CD25-positive, alpha-beta regulatory T cell primary cell |
| cingulate gyrus tissue                                                 |
| duodenal mucosa tissue                                                 |
| Parathyroid adenoma tissue                                             |
| WI38 cell line                                                         |
| angular gyrus tissue                                                   |
| WERI-Rb-1 cell line                                                    |
| adipocyte in vitro differentiated cells                                |
| BE2C cell line                                                         |
| effector memory CD4-positive, alpha-beta T cell primary cell           |
| iPS-18c induced pluripotent stem cell line                             |
| HEK293 cell line                                                       |
| mesodermal cell in vitro differentiated cells                          |
| iPS-15b induced pluripotent stem cell line                             |
| GM12875 cell line                                                      |
| fibroblast of pulmonary artery primary cell                            |
| heart right ventricle tissue                                           |
| iPS-11a induced pluripotent stem cell line                             |
| epithelial cell of esophagus primary cell                              |
| fibroblast of mammary gland primary cell                               |
| AG10803 cell line                                                      |
| choroid plexus epithelial cell primary cell                            |
| AG04449 cell line                                                      |
| aorta tissue                                                           |
| heart tissue                                                           |
| T-helper 17 cell primary cell                                          |
| adipose tissue tissue                                                  |
+------------------------------------------------------------------------+
152 rows in set (34 min 54.65 sec)

+----------+

mysql> select count(*) from tad;
+----------+
| count(*) |
+----------+
|    74424 |
+----------+
1 row in set (0.00 sec)

mysql> select distinct(cell_or_tissue) from tad;
+----------------+
| cell_or_tissue |
+----------------+
| Liver          |
| Thymus         |
| H1-MSC         |
| H1-ESC         |
| H1-MES         |
| H1-NPC         |
| VentricleLeft  |
| H1-TRO         |
| Aorta          |
| K562           |
| HMEC           |
| HUVEC          |
| NHEK           |
| G401           |
| RPMI7951       |
| AdrenalGland   |
| PANC1          |
| LNCaP          |
| Lung           |
| SKNMC          |
| Bladder        |
| Spleen         |
| Cortex         |
| Pancreas       |
| VentricleRight |
| Caki2          |
| Muscle         |
| SKNDZ          |
| Bowel          |
| GM12878        |
| KBM7           |
| IMR90          |
| SKMEL5         |
| A549           |
| NCIH460        |
| T470           |
| SJCRH30        |
+----------------+
37 rows in set (0.41 sec)

+----------+

mysql> select count(*) from tf_binding;
+----------+
| count(*) |
+----------+
| 79812241 |
+----------+
1 row in set (0.00 sec)

mysql> select distinct(tf_name) from tf_binding;
+----------+
| tf_name  |
+----------+
| TAL1     |
| GABPA    |
| FOXP1    |
| NME2     |
| NRF1     |
| CTCF     |
| FOXA1    |
| RAD21    |
| KDM5B    |
| SMC3     |
| SMAD1    |
| STAG1    |
| NKX3-1   |
| AR       |
| E2F6     |
| MAX      |
| HIF1A    |
| MYC      |
| CTCFL    |
| FLI1     |
| ERG      |
| RARA     |
| GATA2    |
| HOXC11   |
| JUND     |
| HNF4A    |
| NIPBL    |
| FOS      |
| TBP      |
| MAFK     |
| HSF1     |
| RUNX1    |
| RELA     |
| NR3C1    |
| VDR      |
| LYL1     |
| MBD2     |
| HDAC2    |
| CDK9     |
| NOTCH1   |
| LMO2     |
| KDM1A    |
| SETDB1   |
| EZH2     |
| SPI1     |
| BRD4     |
| ZBTB33   |
| GTF2B    |
| RXRA     |
| YY1      |
| EBF1     |
| SRF      |
| REST     |
| PBX2     |
| EZH1     |
| PGR      |
| CDK6     |
| CEBPB    |
| ESR1     |
| BHLHE40  |
| CREB1    |
| ZBTB7A   |
| NANOG    |
| NCOR     |
| NR5A2    |
| CHD8     |
| SMARCC1  |
| E2F1     |
| GATA1    |
| EP300    |
| TP53     |
| MED1     |
| RCOR1    |
| CDK2     |
| ZNF750   |
| NFYB     |
| IRF4     |
| RBBP5    |
| TRIM28   |
| KLF4     |
| NR2F2    |
| TAF2     |
| ARNT     |
| SOX2     |
| NCOR2    |
| RELB     |
| GATA6    |
| PPARG    |
| CBFB     |
| CTNNB1   |
| PAX5     |
| PHF8     |
| TOP1     |
| NFATC1   |
| MAF      |
| ZNF263   |
| KDM4C    |
| ZFX      |
| CDK8     |
| PML      |
| TAF1     |
| MYB      |
| CSNK2A1  |
| TAF3     |
| CDK7     |
| SMAD5    |
| RUNX2    |
| SIN3A    |
| AFF4     |
| CREM     |
| KDM4B    |
| TCF12    |
| JMJD1C   |
| HCFC1    |
| CEBPD    |
| HDAC1    |
| SOX9     |
| SP1      |
| POU5F1   |
| NKX2-2   |
| TCF7L2   |
| POU2F2   |
| MXI1     |
| BCL6     |
| TCF3     |
| ZKSCAN1  |
| GFI1B    |
| NEUROD1  |
| TFAP2C   |
| KDM4A    |
| ZHX2     |
| ATF3     |
| ZHX1     |
| NFE2     |
| SIN3B    |
| STAT1    |
| EGLN2    |
| CHD1     |
| NFE2L2   |
| HBP1     |
| CREBBP   |
| ELL2     |
| KLF9     |
| MAFB     |
| THAP1    |
| SNAPC1   |
| ETS1     |
| MAZ      |
| CHD2     |
| JUNB     |
| ASXL1    |
| BRCA1    |
| NFKB1    |
| ZMIZ1    |
| TBL1XR1  |
| ATF1     |
| EGR1     |
| ZNF143   |
| NFYA     |
| GATA3    |
| BRD1     |
| KDM5A    |
| JUN      |
| ZNF175   |
| MAFF     |
| MYCN     |
| RBPJ     |
| UBTF     |
| DNMT3B   |
| CXXC4    |
| STAT5B   |
| SUZ12    |
| ASH2L    |
| SIX5     |
| NKX2-1   |
| CBX1     |
| ELF1     |
| IRF1     |
| E2F4     |
| RFX5     |
| TBL1X    |
| FOXA2    |
| STAT3    |
| GATAD1   |
| KLF5     |
| SP4      |
| RUNX1T1  |
| HMGN3    |
| ZNF639   |
| CCNT2    |
| IRF2     |
| USF1     |
| GTF2F1   |
| DEK      |
| ZNF384   |
| OCA2     |
| FOXK1    |
| TCF7     |
| CBX3     |
| ZBED1    |
| FOSL2    |
| ATF2     |
| BATF     |
| FOXP2    |
| XBP1     |
| ARRB1    |
| SMC1A    |
| RB1      |
| TARDBP   |
| SREBF1   |
| TRIM24   |
| USF2     |
| BCOR     |
| BACH2    |
| ZEB1     |
| TAF7     |
| IKZF1    |
| E2F7     |
| SMAD3    |
| TEAD4    |
| CLOCK    |
| PAX8     |
| NR2C2    |
| ZBTB11   |
| BCLAF1   |
| RUNX3    |
| KLF1     |
| TFAP4    |
| MYOD1    |
| PDX1     |
| NRIP1    |
| NCOA2    |
| BCL3     |
| PBX3     |
| MEF2A    |
| MEF2C    |
| MITF     |
| SNAI2    |
| NR2F1    |
| CBX2     |
| JARID2   |
| PCGF2    |
| TP73     |
| MLLT1    |
| DUX4     |
| MEIS1    |
| TBX21    |
| RNF2     |
| ETV1     |
| BMI1     |
| EED      |
| USP7     |
| PRDM14   |
| SMARCB1  |
| POU4F2   |
| SMARCA4  |
| CBFA2T2  |
| PKNOX1   |
| HEXIM1   |
| KDM6B    |
| TFDP1    |
| RXR      |
| MTA3     |
| MYNN     |
| TP63     |
| MED26    |
| DEAF1    |
| FOSL1    |
| GATA4    |
| FOXO1    |
| ZFP36    |
| ELK1     |
| PIAS1    |
| MBD1     |
| SOX6     |
| NONO     |
| ARID3A   |
| JMJD6    |
| NFIC     |
| E2F5     |
| BCL11A   |
| BACH1    |
| ZNF24    |
| ETV4     |
| TEAD1    |
| NELFA    |
| KLF6     |
| DPF2     |
| MTA2     |
| CTBP2    |
| TWIST1   |
| STAT5A   |
| TRIM22   |
| ID3      |
| SND1     |
| HDGF     |
| BRD2     |
| BRD3     |
| DAXX     |
| ATF7     |
| ZNF207   |
| AHR      |
| PROX1    |
| GRHL2    |
| MEIS2    |
| ZEB2     |
| CUX1     |
| KLF13    |
| SREBF2   |
| BHLHE22  |
| BMPR1A   |
| ASCL1    |
| ZNF92    |
| ELF3     |
| HNF1B    |
| MNT      |
| MAML3    |
| FOXM1    |
| KLF3     |
| SIX2     |
| ESRRA    |
| RYBP     |
| HOXB13   |
| MBD3     |
| HDAC8    |
| NCOA3    |
| CTBP1    |
| HOXA9    |
| NCOR1    |
| TFAP2A   |
| NR2F6    |
| SMAD2    |
| SMAD4    |
| NR1H2    |
| HNF4G    |
| CDX2     |
| NR4A1    |
| GTF3C2   |
| ZMYND11  |
| MTA1     |
| NSD2     |
| SFMBT1   |
| ARID1A   |
| DDX5     |
| SAP30    |
| ONECUT1  |
| EOMES    |
| TCF4     |
| T        |
| ESR2     |
| MIER1    |
| HINFP    |
| NBN      |
| ZNF740   |
| CHD7     |
| ARID1B   |
| TLE3     |
| LHX2     |
| UBN1     |
| SIRT6    |
| SSRP1    |
| SPIB     |
| CRY1     |
| L3MBTL2  |
| ARNTL    |
| RBP2     |
| FOXK2    |
| KMT2B    |
| ELK4     |
| SP2      |
| ZC3H11A  |
| FANCL    |
| CREB3L1  |
| SOX11    |
| NCOA1    |
| NFRKB    |
| PRKDC    |
| NFKB2    |
| ETV6     |
| ARID2    |
| TGIF2    |
| CHD4     |
| CEBPZ    |
| PALB2    |
| ZNF274   |
| DIDO1    |
| CREB3    |
| CEBPG    |
| ZC3H8    |
| ZNF584   |
| ZNF589   |
| CEBPA    |
| RAD51    |
| ZNF318   |
| STAT2    |
| TCF25    |
| ADNP     |
| ZNF644   |
| IRF9     |
| SOX13    |
| MAFG     |
| HMBOX1   |
| LDB1     |
| MEF2B    |
| POU2F1   |
| HDAC6    |
| YAP1     |
| TERF1    |
| IRF3     |
| RFX1     |
| KAT2B    |
| ZNF217   |
| EPAS1    |
| RFX2     |
| ZNF197   |
| SMARCC2  |
| ZNF83    |
| ATRX     |
| KMT2A    |
| DDX20    |
| TSC22D4  |
| HHEX     |
| PTTG1    |
| PTRF     |
| PYGO2    |
| TERF2    |
| GSPT2    |
| RUNX     |
| MYBL2    |
| ZBTB16   |
| THAP11   |
| PAX6     |
| BRF2     |
| NFKBIA   |
| SATB1    |
| SOX10    |
| TEAD2    |
| NFE2L1   |
| ILK      |
| CIITA    |
| OTX2     |
| HDAC3    |
| NR1H3    |
| LEF1     |
| ZNF766   |
| PPARGC1A |
| ZNF165   |
| ICE1     |
| ICE2     |
| BAHD1    |
| SNAPC4   |
| FOXJ2    |
| ATF4     |
| MBD4     |
| GLYR1    |
| PRDM1    |
| CBX4     |
| HES1     |
| CARM1    |
| HOXB7    |
| SVIL     |
| ELL      |
| SMC4     |
| BDP1     |
| KAT2A    |
| C17ORF49 |
| TRRAP    |
| SGF29    |
| MLLT3    |
| ELF5     |
| BRF1     |
| EHF      |
| ZNF592   |
| EWSR1    |
| ZZZ3     |
| SIRT3    |
| WRNIP1   |
| SUPT20H  |
| ZNF335   |
| MDM2     |
| OVOL2    |
| LRWD1    |
| KDM1B    |
| DPPA3    |
| TET3     |
+----------+
485 rows in set (5 min 0.13 sec)
```
