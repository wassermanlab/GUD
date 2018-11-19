#!/usr/bin/env python

import os, sys, re
import argparse
from binning import assign_bin
from datetime import date
import getpass
import pybedtools
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
from urllib2 import unquote
import warnings

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.tss import TSS

#-------------#
# Definitions #
#-------------#

sample_names = {
    "CNhs11844": "HCC1806",
    "CNhs11251": "BALL-1",
    "CNhs11282": "NALM-6",
    "CNhs10746": "HPB-ALL",
    "CNhs11253": "Jurkat",
    "CNhs13053": "KG-1",
    "CNhs13054": "HYT-1",
    "CNhs13502": "Kasumi-1",
    "CNhs13052": "Kasumi-6",
    "CNhs11864": "NKM-1",
    "CNhs13055": "HL60",
    "CNhs13503": "FKH-1",
    "CNhs13504": "HNT-34",
    "CNhs13056": "EoL-1",
    "CNhs13057": "EoL-3",
    "CNhs13050": "NOMO-1",
    "CNhs13051": "P31/FUJ",
    "CNhs10722": "THP-1",
    "CNhs13058": "U-937 DE-4",
    "CNhs13059": "EEB",
    "CNhs13060": "F-36E",
    "CNhs13505": "F-36P",
    "CNhs11888": "MKPL-1",
    "CNhs13049": "M-MOK",
    "CNhs11882": "IM95m",
    "CNhs11051": "adipocyte (breast)",
    "CNhs11969": "adipocyte (breast)",
    "CNhs11054": "adipocyte (omentum)",
    "CNhs12067": "adipocyte (omentum)",
    "CNhs12068": "adipocyte (omentum)",
    "CNhs12069": "adipocyte (perirenal)",
    "CNhs12494": "adipocyte (subcutaneous)",
    "CNhs11371": "adipocyte (subcutaneous)",
    "CNhs12017": "adipocyte (subcutaneous)",
    "CNhs10615": "adipose tissue",
    "CNhs11893": "SW-13",
    "CNhs10738": "ATN-1",
    "CNhs11838": "SW 1573",
    "CNhs11325": "epithelial cells (alveolus)",
    "CNhs12084": "epithelial cells (alveolus)",
    "CNhs11341": "epithelial cells (amniotic membrane)",
    "CNhs12125": "epithelial cells (amniotic membrane)",
    "CNhs12502": "amniotic membrane cells",
    "CNhs12503": "amniotic membrane cells",
    "CNhs12379": "amniotic membrane cells",
    "CNhs13793": "amygdala",
    "CNhs12311": "amygdala",
    "CNhs10745": "8305C",
    "CNhs11881": "Ki-JK",
    "CNhs11889": "RPMI 2650",
    "CNhs10876": "anulus pulposus cells",
    "CNhs12064": "anulus pulposus cells",
    "CNhs11760": "aorta",
    "CNhs12842": "appendix",
    "CNhs11725": "TC-YIK",
    "CNhs11321": "astrocytes (cerebellum)",
    "CNhs12081": "astrocytes (cerebellum)",
    "CNhs12117": "astrocytes (cerebellum)",
    "CNhs10864": "astrocytes (cerebral cortex)",
    "CNhs11960": "astrocytes (cerebral cortex)",
    "CNhs12005": "astrocytes (cerebral cortex)",
    "CNhs10742": "TM-31",
    "CNhs11932": "TE 354.T",
    "CNhs12575": "basophils",
    "CNhs10744": "RPMI1788",
    "CNhs10750": "HuCCT1",
    "CNhs11265": "TFK-1",
    "CNhs11845": "MV-4-11",
    "CNhs10616": "bladder",
    "CNhs11761": "blood",
    "CNhs12331": "GM12878",
    "CNhs12332": "GM12878",
    "CNhs12333": "GM12878",
    "CNhs11931": "StromaNKtert",
    "CNhs11796": "brain",
    "CNhs10617": "brain",
    "CNhs11797": "brain (fetal)",
    "CNhs11792": "breast",
    "CNhs11943": "MCF-7",
    "CNhs10736": "MDA-MB-453",
    "CNhs12054": "epithelial cells (bronchi)",
    "CNhs12058": "epithelial cells (bronchi)",
    "CNhs12062": "epithelial cells (bronchi)",
    "CNhs11862": "KNS-62",
    "CNhs11840": "NCI-H358",
    "CNhs11841": "ChaGo-K-1",
    "CNhs10739": "DAUDI",
    "CNhs11268": "RAJI",
    "CNhs11834": "NCI-H1770",
    "CNhs11846": "SK-PN-DW",
    "CNhs11747": "JHUCS-1",
    "CNhs12341": "myocytes (heart)",
    "CNhs12350": "myocytes (heart)",
    "CNhs12571": "myocytes (heart)",
    "CNhs13802": "caudate nucleus",
    "CNhs12321": "caudate nucleus",
    "CNhs13224": "monocytes (CD14-positive CD16-negative)",
    "CNhs13541": "monocytes (CD14-positive CD16-positive)",
    "CNhs13207": "monocytes (CD14-negative CD16-positive)",
    "CNhs13208": "monocytes (CD14-positive CD16-positive)",
    "CNhs13216": "monocytes (CD14-positive CD16-negative)",
    "CNhs13540": "monocytes (CD14-positive CD16-negative)",
    "CNhs13548": "monocytes (CD14-negative CD16-positive)",
    "CNhs13549": "monocytes (CD14-positive CD16-positive)",
    "CNhs10858": "endothelial progenitor cells (derived from CD14-positive monocytes)",
    "CNhs11897": "endothelial progenitor cells (derived from CD14-positive monocytes)",
    "CNhs11904": "endothelial progenitor cells (derived from CD14-positive monocytes)",
    "CNhs10852": "monocytes (CD14-positive)",
    "CNhs11954": "monocytes (CD14-positive)",
    "CNhs11997": "monocytes (CD14-positive)",
    "CNhs13468": "monocytes (CD14-positive)",
    "CNhs13484": "monocytes (CD14-positive)",
    "CNhs13491": "monocytes (CD14-positive)",
    "CNhs12343": "B cells (CD19-positive)",
    "CNhs12352": "B cells (CD19-positive)",
    "CNhs12354": "B cells (CD19-positive)",
    "CNhs12588": "CD34-positive stem cells (derived from bone marrow)",
    "CNhs13539": "T cells (CD4-positive CD25-negative CD45RA-negative memory conventional)",
    "CNhs13195": "T cells (CD4-positive CD25-positive CD45RA-negative memory regulatory)",
    "CNhs13206": "T cells (CD4-positive CD25-positive CD45RA-negative memory regulatory)",
    "CNhs13538": "T cells (CD4-positive CD25-positive CD45RA-negative memory regulatory)",
    "CNhs13223": "T cells (CD4-positive CD25-negative CD45RA-positive naive conventional)",
    "CNhs13205": "T cells (CD4-positive CD25-negative CD45RA-positive naive conventional)",
    "CNhs13512": "T cells (CD4-positive CD25-negative CD45RA-positive naive conventional)",
    "CNhs13513": "T cells (CD4-positive CD25-positive CD45RA-positive naive regulatory)",
    "CNhs10853": "T cells (CD4-positive)",
    "CNhs11955": "T cells (CD4-positive)",
    "CNhs11998": "T cells (CD4-positive)",
    "CNhs10854": "T cells (CD8-positive)",
    "CNhs11956": "T cells (CD8-positive)",
    "CNhs11999": "T cells (CD8-positive)",
    "CNhs13799": "cerebellum",
    "CNhs12323": "cerebellum",
    "CNhs11795": "cerebellum",
    "CNhs12840": "meninges (brain)",
    "CNhs11288": "D98-AH2",
    "CNhs11289": "ME-180",
    "CNhs10618": "cervix",
    "CNhs11283": "HuH-28",
    "CNhs11923": "chondrocytes (de-differentiated)",
    "CNhs11372": "chondrocytes (de-differentiated)",
    "CNhs12020": "chondrocytes (de-differentiated)",
    "CNhs11373": "chondrocytes (re-differentiated)",
    "CNhs12021": "chondrocytes (re-differentiated)",
    "CNhs10740": "BeWo",
    "CNhs11875": "SCH",
    "CNhs11820": "T3M-3",
    "CNhs12504": "chorionic membrane cells",
    "CNhs12506": "chorionic membrane cells",
    "CNhs12380": "chorionic membrane cells",
    "CNhs11714": "SKW-3",
    "CNhs11886": "KCL-22",
    "CNhs11250": "K562",
    "CNhs12334": "K562",
    "CNhs12335": "K562",
    "CNhs12336": "K562",
    "CNhs10727": "KU812",
    "CNhs11865": "MEG-A2",
    "CNhs10871": "ciliary epithelial cells",
    "CNhs11966": "ciliary epithelial cells",
    "CNhs12009": "ciliary epithelial cells",
    "CNhs11745": "JHOC-5",
    "CNhs11930": "TEN",
    "CNhs11794": "colon",
    "CNhs10619": "colon",
    "CNhs11280": "CACO-2",
    "CNhs10737": "COLO-320",
    "CNhs11780": "colon (fetal)",
    "CNhs11045": "COBL-a",
    "CNhs11336": "epithelial cells (cornea)",
    "CNhs12123": "epithelial cells (cornea)",
    "CNhs10649": "corpus callosum",
    "CNhs10855": "dendritic cells (derived from immature monocytes)",
    "CNhs11062": "dendritic cells (derived from immature monocytes)",
    "CNhs12000": "dendritic cells (derived from immature monocytes)",
    "CNhs10857": "plasmacytoid dendritic cells",
    "CNhs11779": "diaphragm (fetal)",
    "CNhs12610": "diencephalon",
    "CNhs11741": "CTB-1",
    "CNhs11100": "KLM-1",
    "CNhs11259": "MIA Paca2",
    "CNhs12846": "ductus deferens",
    "CNhs11781": "duodenum (fetal)",
    "CNhs12996": "duodenum (fetal)",
    "CNhs10648": "dura mater",
    "CNhs11046": "HEK293/SLAM",
    "CNhs11731": "1B2C6",
    "CNhs11732": "1C3D3",
    "CNhs11733": "1C3IKEI",
    "CNhs11814": "2C6",
    "CNhs11266": "OMC-2",
    "CNhs11249": "OMC-9",
    "CNhs11748": "JHUEM-1",
    "CNhs10837": "endothelial cells (aorta)",
    "CNhs12495": "endothelial cells (aorta)",
    "CNhs11375": "endothelial cells (aorta)",
    "CNhs12022": "endothelial cells (aorta)",
    "CNhs12496": "endothelial cells (artery)",
    "CNhs11977": "endothelial cells (artery)",
    "CNhs12023": "endothelial cells (artery)",
    "CNhs10865": "endothelial cells (lymph node)",
    "CNhs11901": "endothelial cells (lymph node)",
    "CNhs11906": "endothelial cells (lymph node)",
    "CNhs11925": "endothelial cells (microvasculature)",
    "CNhs11376": "endothelial cells (microvasculature)",
    "CNhs12024": "endothelial cells (microvasculature)",
    "CNhs11926": "endothelial cells (thorax)",
    "CNhs11978": "endothelial cells (thorax)",
    "CNhs10872": "endothelial cells (umbilical vein)",
    "CNhs11967": "endothelial cells (umbilical vein)",
    "CNhs12010": "endothelial cells (umbilical vein)",
    "CNhs12497": "endothelial cells (vein)",
    "CNhs11377": "endothelial cells (vein)",
    "CNhs12026": "endothelial cells (vein)",
    "CNhs10743": "A431",
    "CNhs10748": "Ca Ski",
    "CNhs12847": "epididymis",
    "CNhs11247": "HS-ES-1",
    "CNhs12325": "HeLa-S3",
    "CNhs12326": "HeLa-S3",
    "CNhs12327": "HeLa-S3",
    "CNhs11323": "epithelial cells (esophagus)",
    "CNhs10620": "esophagus",
    "CNhs11836": "Hs 863.T",
    "CNhs10728": "H-EMC-SS",
    "CNhs11762": "eye (fetal)",
    "CNhs10874": "fibroblast (aortic adventitia)",
    "CNhs12011": "fibroblast (aortic adventitia)",
    "CNhs12498": "fibroblast (heart)",
    "CNhs11378": "fibroblast (heart)",
    "CNhs12027": "fibroblast (heart)",
    "CNhs11909": "fibroblast (heart)",
    "CNhs12057": "fibroblast (heart)",
    "CNhs12061": "fibroblast (heart)",
    "CNhs11319": "fibroblast (choroid plexus)",
    "CNhs12344": "fibroblast (choroid plexus)",
    "CNhs11339": "fibroblast (conjunctiva)",
    "CNhs12499": "fibroblast (dermis)",
    "CNhs11379": "fibroblast (dermis)",
    "CNhs12028": "fibroblast (dermis)",
    "CNhs12052": "fibroblast (dermis)",
    "CNhs12055": "fibroblast (dermis)",
    "CNhs12059": "fibroblast (dermis)",
    "CNhs10866": "fibroblast (gingiva)",
    "CNhs11961": "fibroblast (gingiva)",
    "CNhs12006": "fibroblast (gingiva)",
    "CNhs10848": "fibroblast (gingiva)",
    "CNhs11952": "fibroblast (gingiva)",
    "CNhs11322": "fibroblast (lymph node)",
    "CNhs12118": "fibroblast (lymph node)",
    "CNhs10867": "fibroblast (periodontal ligament)",
    "CNhs11962": "fibroblast (periodontal ligament)",
    "CNhs11907": "fibroblast (periodontal ligament)",
    "CNhs12493": "fibroblast (periodontal ligament)",
    "CNhs11953": "fibroblast (periodontal ligament)",
    "CNhs11996": "fibroblast (periodontal ligament)",
    "CNhs10878": "fibroblast (pulmonary artery)",
    "CNhs11351": "fibroblast (skin)",
    "CNhs11914": "fibroblast (skin)",
    "CNhs11860": "HT-1080",
    "CNhs11842": "GCT TIB-223",
    "CNhs10647": "frontal lobe",
    "CNhs12848": "gall bladder",
    "CNhs11256": "TGBC14TKB",
    "CNhs10733": "TGBC2TKB",
    "CNhs11737": "MKN1",
    "CNhs11819": "MKN45",
    "CNhs11286": "AZ521",
    "CNhs11738": "ECC12",
    "CNhs11274": "LU65",
    "CNhs10751": "Lu99B",
    "CNhs11061": "epithelial cells (gingiva)",
    "CNhs11896": "epithelial cells (gingiva)",
    "CNhs11903": "epithelial cells (gingiva)",
    "CNhs11824": "HOKUG",
    "CNhs11185": "A172",
    "CNhs11248": "A172",
    "CNhs11272": "T98G",
    "CNhs10731": "GI-1",
    "CNhs13801": "globus pallidus",
    "CNhs12319": "globus pallidus",
    "CNhs11740": "KGN",
    "CNhs12501": "hair follicle dermal papilla cells",
    "CNhs11979": "hair follicle dermal papilla cells",
    "CNhs12030": "hair follicle dermal papilla cells",
    "CNhs12339": "hair follicle outer root sheath cells",
    "CNhs12347": "hair follicle outer root sheath cells",
    "CNhs11843": "Mo",
    "CNhs10621": "heart",
    "CNhs10653": "heart (fetal)",
    "CNhs12855": "mitral valve",
    "CNhs12856": "pulmonic valve",
    "CNhs12857": "tricuspid valve",
    "CNhs13479": "HEp-2",
    "CNhs13500": "HEp-2",
    "CNhs13501": "HEp-2",
    "CNhs11868": "LI90",
    "CNhs12075": "endothelial cells (liver sinusoid)",
    "CNhs12092": "endothelial cells (liver sinusoid)",
    "CNhs11335": "hepatic stellate cells",
    "CNhs12093": "hepatic stellate cells",
    "CNhs11742": "HuH-6",
    "CNhs12328": "HepG2",
    "CNhs12329": "HepG2",
    "CNhs12330": "HepG2",
    "CNhs12340": "hepatocytes",
    "CNhs12349": "hepatocytes",
    "CNhs12626": "hepatocytes",
    "CNhs11271": "Li-7",
    "CNhs11891": "WIL2-NS",
    "CNhs13795": "hippocampus",
    "CNhs12312": "hippocampus",
    "CNhs11715": "HD-Mar2",
    "CNhs13537": "langerhans cells (immature)",
    "CNhs13480": "langerhans cells (immature)",
    "CNhs10646": "insula",
    "CNhs10875": "epithelial cells (intestine, polarized)",
    "CNhs12596": "iris pigment epithelium",
    "CNhs11064": "keratinocytes (epidermis)",
    "CNhs11381": "keratinocytes (epidermis)",
    "CNhs12031": "keratinocytes (epidermis)",
    "CNhs10879": "keratinocytes (oral)",
    "CNhs11880": "HKA-1",
    "CNhs11337": "keratocytes",
    "CNhs12095": "keratocytes",
    "CNhs10622": "kidney",
    "CNhs10652": "kidney (fetal)",
    "CNhs11277": "IA-LM",
    "CNhs12806": "NCI-H460",
    "CNhs11825": "SKG-II-SF",
    "CNhs11790": "left atrium",
    "CNhs11789": "left ventricle",
    "CNhs11848": "G-402",
    "CNhs11722": "10964C",
    "CNhs11723": "15242A",
    "CNhs11724": "15425",
    "CNhs11750": "SRA 01/04",
    "CNhs12342": "epithelial cells (lens)",
    "CNhs12568": "epithelial cells (lens)",
    "CNhs12572": "epithelial cells (lens)",
    "CNhs11859": "MEG-01",
    "CNhs11870": "KMLS-1",
    "CNhs11851": "SW 872",
    "CNhs10624": "liver",
    "CNhs11798": "liver (fetal)",
    "CNhs13808": "locus coeruleus",
    "CNhs12322": "locus coeruleus",
    "CNhs11275": "A549",
    "CNhs10726": "PC-14",
    "CNhs10625": "lung",
    "CNhs11680": "lung (fetal)",
    "CNhs11786": "lung (right lower lobe)",
    "CNhs11852": "DS-1",
    "CNhs11788": "lymph node",
    "CNhs11935": "MLMA",
    "CNhs10861": "macrophages (derived from monocytes)",
    "CNhs11899": "macrophages (derived from monocytes)",
    "CNhs12003": "macrophages (derived from monocytes)",
    "CNhs10730": "DJM-1",
    "CNhs13550": "Malassez epithelial cells",
    "CNhs13551": "Malassez epithelial cells",
    "CNhs11077": "epithelial cells (breast)",
    "CNhs11382": "epithelial cells (breast)",
    "CNhs12032": "epithelial cells (breast)",
    "CNhs12566": "mast cells",
    "CNhs12594": "mast cells",
    "CNhs12593": "mast cells",
    "CNhs12592": "mast cells",
    "CNhs10732": "HSQ-89",
    "CNhs13796": "medial frontal gyrus",
    "CNhs13809": "medial temporal gyrus",
    "CNhs12310": "medial temporal gyrus",
    "CNhs13800": "medulla oblongata",
    "CNhs12315": "medulla oblongata",
    "CNhs10645": "medulla oblongata",
    "CNhs12805": "D283 Med",
    "CNhs11861": "ONS-76",
    "CNhs12570": "melanocytes (dark donor)",
    "CNhs11303": "melanocytes (light donor)",
    "CNhs11383": "melanocytes (light donor)",
    "CNhs12033": "melanocytes (light donor)",
    "CNhs11281": "COLO 679",
    "CNhs11254": "G-361",
    "CNhs11320": "meningeal cells",
    "CNhs12080": "meningeal cells",
    "CNhs12731": "meningeal cells",
    "CNhs11945": "HKBMM",
    "CNhs12838": "MKL-1",
    "CNhs12839": "MS-1",
    "CNhs12363": "mesenchymal precursor cells (adipose tissue)",
    "CNhs12364": "mesenchymal precursor cells (adipose tissue)",
    "CNhs12365": "mesenchymal precursor cells (adipose tissue)",
    "CNhs12366": "mesenchymal precursor cells (bone marrow)",
    "CNhs12367": "mesenchymal precursor cells (bone marrow)",
    "CNhs13098": "mesenchymal precursor cells (bone marrow)",
    "CNhs12368": "mesenchymal precursor cells (heart)",
    "CNhs12369": "mesenchymal precursor cells (heart)",
    "CNhs12370": "mesenchymal precursor cells (heart)",
    "CNhs12371": "mesenchymal precursor cells (heart)",
    "CNhs11718": "Hu5/E18",
    "CNhs10844": "mesenchymal stem cells (adipose tissue)",
    "CNhs11345": "mesenchymal stem cells (adipose tissue)",
    "CNhs12922": "mesenchymal stem cells (adipose tissue)",
    "CNhs11349": "mesenchymal stem cells (amniotic membrane)",
    "CNhs12104": "mesenchymal stem cells (amniotic membrane)",
    "CNhs11344": "mesenchymal stem cells (bone marrow)",
    "CNhs12100": "mesenchymal stem cells (bone marrow)",
    "CNhs12126": "mesenchymal stem cells (bone marrow)",
    "CNhs10845": "mesenchymal stem cells (liver)",
    "CNhs12730": "mesenchymal stem cells (liver)",
    "CNhs12492": "mesenchymal stem cells (umbilical cord)",
    "CNhs11347": "mesenchymal stem cells (umbilical cord)",
    "CNhs12127": "mesenchymal stem cells (umbilical cord)",
    "CNhs10850": "mesothelial cells",
    "CNhs12012": "mesothelial cells",
    "CNhs11263": "ACC-MESO-1",
    "CNhs11264": "ACC-MESO-4",
    "CNhs13066": "Mero-25",
    "CNhs13067": "Mero-41",
    "CNhs13068": "Mero-48a",
    "CNhs13069": "Mero-82",
    "CNhs13070": "Mero-83",
    "CNhs13072": "Mero-84",
    "CNhs13073": "Mero-95",
    "CNhs13063": "NCI-H2052",
    "CNhs13062": "NCI-H226",
    "CNhs13064": "NCI-H2452",
    "CNhs13061": "NCI-H28",
    "CNhs13074": "No36",
    "CNhs13075": "ONE58",
    "CNhs12316": "middle temporal gyrus",
    "CNhs13535": "langerhans cells (migratory)",
    "CNhs13536": "langerhans cells (migratory)",
    "CNhs13547": "langerhans cells (migratory)",
    "CNhs11944": "HTMMT",
    "CNhs11752": "JHOM-1",
    "CNhs11873": "MCAS",
    "CNhs11350": "multipotent stem cells (umbilical cord blood)",
    "CNhs12105": "multipotent stem cells (umbilical cord blood)",
    "CNhs11858": "HuT 102 TIB-162",
    "CNhs11934": "SKM-1",
    "CNhs11258": "PCM6",
    "CNhs10870": "myoblasts",
    "CNhs11965": "myoblasts",
    "CNhs11908": "myoblasts",
    "CNhs11729": "MFH-ino",
    "CNhs11821": "NMFH-1",
    "CNhs12589": "epithelial cells (nasal)",
    "CNhs12574": "epithelial cells (nasal)",
    "CNhs10859": "natural killer cells",
    "CNhs11957": "natural killer cells",
    "CNhs12001": "natural killer cells",
    "CNhs11063": "neural stem cells",
    "CNhs11384": "neural stem cells",
    "CNhs11276": "CHP-134",
    "CNhs11284": "NB-1",
    "CNhs11818": "NBsusSR",
    "CNhs11811": "NH-12",
    "CNhs11744": "FU-RPNT-1",
    "CNhs11753": "FU-RPNT-2",
    "CNhs11866": "TASK1",
    "CNhs11853": "SK-N-MC",
    "CNhs11854": "Hs 53.T",
    "CNhs12338": "neurons",
    "CNhs12726": "neurons",
    "CNhs13815": "neurons",
    "CNhs10862": "neutrophils",
    "CNhs11959": "neutrophils",
    "CNhs11905": "neutrophils",
    "CNhs11867": "KHYG-1",
    "CNhs10747": "P30/OHK",
    "CNhs11894": "HEPM",
    "CNhs11950": "FHs 74 Int",
    "CNhs10644": "nucleus accumbens",
    "CNhs10881": "nucleus pulposus cells",
    "CNhs12019": "nucleus pulposus cells",
    "CNhs12063": "nucleus pulposus cells",
    "CNhs13798": "occipital cortex",
    "CNhs12320": "occipital cortex",
    "CNhs11787": "occipital lobe",
    "CNhs11784": "occipital lobe (fetal)",
    "CNhs10643": "occipital pole",
    "CNhs13816": "olfactory epithelium",
    "CNhs13817": "olfactory epithelium",
    "CNhs13818": "olfactory epithelium",
    "CNhs13819": "olfactory epithelium",
    "CNhs12611": "olfactory region",
    "CNhs13449": "optic nerve",
    "CNhs10752": "Ca9-22",
    "CNhs11287": "HO-1-u-1",
    "CNhs11717": "HSC-3",
    "CNhs11810": "SAS",
    "CNhs11311": "osteoblasts (differentiated)",
    "CNhs11980": "osteoblasts (differentiated)",
    "CNhs12035": "osteoblasts (differentiated)",
    "CNhs11385": "osteoblasts",
    "CNhs12036": "osteoblasts",
    "CNhs11835": "Hs 706.T",
    "CNhs11279": "143B/TK^(-)neo^(R)",
    "CNhs11290": "HS-Os-1",
    "CNhs10626": "ovary",
    "CNhs11856": "Hs 925.T",
    "CNhs11756": "pancreas",
    "CNhs11832": "NOR-P1",
    "CNhs10877": "stromal cells (pancreas)",
    "CNhs11716": "8505C",
    "CNhs10734": "TGBC18TKB",
    "CNhs10642": "paracentral gyrus",
    "CNhs13797": "parietal lobe",
    "CNhs12317": "parietal lobe",
    "CNhs10641": "parietal lobe",
    "CNhs11782": "parietal lobe (fetal)",
    "CNhs12849": "parotid gland",
    "CNhs12850": "penis",
    "CNhs11317": "pericytes",
    "CNhs12079": "pericytes",
    "CNhs10860": "peripheral blood mononuclear cells",
    "CNhs11958": "peripheral blood mononuclear cells",
    "CNhs12002": "peripheral blood mononuclear cells",
    "CNhs11830": "KU-SN",
    "CNhs11849": "Detroit 562",
    "CNhs13804": "pineal gland",
    "CNhs12228": "pineal gland",
    "CNhs13805": "pituitary gland",
    "CNhs12229": "pituitary gland",
    "CNhs10627": "placenta",
    "CNhs11079": "epithelial cells (placenta)",
    "CNhs11386": "epithelial cells (placenta)",
    "CNhs12037": "epithelial cells (placenta)",
    "CNhs12807": "ARH-77",
    "CNhs11933": "SNU-387",
    "CNhs10640": "pons",
    "CNhs10638": "postcentral gyrus",
    "CNhs11052": "preadipocytes (breast)",
    "CNhs11971": "preadipocytes (breast)",
    "CNhs11065": "preadipocytes (omentum)",
    "CNhs11902": "preadipocytes (omentum)",
    "CNhs12013": "preadipocytes (omentum)",
    "CNhs12065": "preadipocytes (perirenal)",
    "CNhs11981": "preadipocytes (subcutaneous)",
    "CNhs12038": "preadipocytes (subcutaneous)",
    "CNhs11082": "preadipocytes (viscera)",
    "CNhs11982": "preadipocytes (viscera)",
    "CNhs12039": "preadipocytes (viscera)",
    "CNhs10628": "prostate",
    "CNhs11260": "DU145",
    "CNhs11243": "PC-3",
    "CNhs11972": "epithelial cells (prostate)",
    "CNhs12014": "epithelial cells (prostate)",
    "CNhs10882": "epithelial cells (prostate, polarized)",
    "CNhs10883": "stromal cells (prostate)",
    "CNhs11973": "stromal cells (prostate)",
    "CNhs12015": "stromal cells (prostate)",
    "CNhs12324": "putamen",
    "CNhs11255": "TT1TKB",
    "CNhs11777": "rectum (fetal)",
    "CNhs10729": "OS-RC-2",
    "CNhs11257": "TUHR10TKB",
    "CNhs11331": "epithelial cells (renal cortex)",
    "CNhs12728": "epithelial cells (renal cortex)",
    "CNhs11332": "epithelial cells (renal)",
    "CNhs12088": "epithelial cells (renal)",
    "CNhs12732": "epithelial cells (renal)",
    "CNhs12074": "endothelial cells (renal glomerulus)",
    "CNhs12086": "endothelial cells (renal glomerulus)",
    "CNhs12624": "endothelial cells (renal glomerulus)",
    "CNhs13080": "endothelial cells (renal glomerulus)",
    "CNhs11333": "mesangial cells",
    "CNhs12121": "mesangial cells",
    "CNhs11330": "epithelial cell (renal proximal tubule)",
    "CNhs12087": "epithelial cell (renal proximal tubule)",
    "CNhs12120": "epithelial cell (renal proximal tubule)",
    "CNhs13552": "reticulocytes",
    "CNhs13553": "reticulocytes",
    "CNhs10636": "retina",
    "CNhs10842": "retinal pigment epithelium",
    "CNhs11338": "retinal pigment epithelium",
    "CNhs12733": "retinal pigment epithelium",
    "CNhs11267": "Y79",
    "CNhs11877": "KYM-1",
    "CNhs11269": "RMS-YM",
    "CNhs11829": "HTST",
    "CNhs12810": "acinar cells (salivary)",
    "CNhs12811": "acinar cells (salivary)",
    "CNhs12812": "acinar cells (salivary)",
    "CNhs11677": "salivary gland",
    "CNhs11183": "HS-PSS",
    "CNhs11245": "HS-PSS",
    "CNhs10847": "sebocytes",
    "CNhs11951": "sebocytes",
    "CNhs12851": "seminal vesicle",
    "CNhs11746": "JHOS-2",
    "CNhs13099": "SK-OV-3-R",
    "CNhs11827": "HTOA",
    "CNhs10851": "Sertoli cells",
    "CNhs10753": "Kato III",
    "CNhs11270": "NUGC-4",
    "CNhs10629": "skeletal muscle",
    "CNhs11084": "myotubes (differentiated from skeletal muscle cells)",
    "CNhs11083": "skeletal muscle cells",
    "CNhs12053": "skeletal muscle cells",
    "CNhs12056": "skeletal muscle cells",
    "CNhs12060": "skeletal muscle cells",
    "CNhs11776": "skeletal muscle (fetal)",
    "CNhs10869": "skeletal muscle satellite cells",
    "CNhs11964": "skeletal muscle satellite cells",
    "CNhs12008": "skeletal muscle satellite cells",
    "CNhs13454": "skeletal muscle (soleus)",
    "CNhs11774": "skin (fetal)",
    "CNhs10884": "epithelial cells (small airways)",
    "CNhs11975": "epithelial cells (small airways)",
    "CNhs12016": "epithelial cells (small airways)",
    "CNhs11885": "HCSC-1",
    "CNhs11736": "ECC10",
    "CNhs11734": "ECC4",
    "CNhs12808": "DMS 144",
    "CNhs11285": "LK-2",
    "CNhs12809": "NCI-H82",
    "CNhs11812": "WA-hT",
    "CNhs10630": "small intestine",
    "CNhs11773": "small intestine (fetal)",
    "CNhs11755": "smooth muscle",
    "CNhs10838": "smooth muscle cells (aorta)",
    "CNhs11085": "smooth muscle cells (aorta)",
    "CNhs11305": "smooth muscle cells (aorta)",
    "CNhs11309": "smooth muscle cells (aorta)",
    "CNhs11086": "smooth muscle cells (brachiocephalic artery)",
    "CNhs12043": "smooth muscle cells (brachiocephalic artery)",
    "CNhs10863": "smooth muscle cells (brain vasculature)",
    "CNhs11900": "smooth muscle cells (brain vasculature)",
    "CNhs12004": "smooth muscle cells (brain vasculature)",
    "CNhs11328": "smooth muscle cells (bronchi)",
    "CNhs12348": "smooth muscle cells (bronchi)",
    "CNhs11087": "smooth muscle cells (carotid)",
    "CNhs12044": "smooth muscle cells (carotid)",
    "CNhs10868": "smooth muscle cells (colon)",
    "CNhs11963": "smooth muscle cells (colon)",
    "CNhs12007": "smooth muscle cells (colon)",
    "CNhs11088": "smooth muscle cells (coronary artery)",
    "CNhs11987": "smooth muscle cells (coronary artery)",
    "CNhs12045": "smooth muscle cells (coronary artery)",
    "CNhs11324": "smooth muscle cells (esophagus)",
    "CNhs11988": "smooth muscle cells (internal thoracic artery)",
    "CNhs12046": "smooth muscle cells (internal thoracic artery)",
    "CNhs11920": "smooth muscle cells (prostate)",
    "CNhs11976": "smooth muscle cells (prostate)",
    "CNhs11989": "smooth muscle cells (pulmonary artery)",
    "CNhs11090": "smooth muscle cells (subclavian artery)",
    "CNhs11990": "smooth muscle cells (subclavian artery)",
    "CNhs12048": "smooth muscle cells (subclavian artery)",
    "CNhs11329": "smooth muscle cells (trachea)",
    "CNhs12894": "smooth muscle cells (trachea)",
    "CNhs10839": "smooth muscle cells (umbilical artery)",
    "CNhs11091": "smooth muscle cells (umbilical artery)",
    "CNhs11991": "smooth muscle cells (umbilical artery)",
    "CNhs12049": "smooth muscle cells (umbilical artery)",
    "CNhs12597": "smooth muscle cells (umbilical vein)",
    "CNhs12569": "smooth muscle cells (umbilical vein)",
    "CNhs11927": "smooth muscle cells (uterus)",
    "CNhs11869": "QGP-1",
    "CNhs13807": "spinal cord",
    "CNhs12227": "spinal cord",
    "CNhs11764": "spinal cord (fetal)",
    "CNhs11857": "Hs 132.T",
    "CNhs10631": "spleen",
    "CNhs10651": "spleen (fetal)",
    "CNhs10741": "SLVL",
    "CNhs11252": "EC-GI-10",
    "CNhs11739": "T3M-5",
    "CNhs11273": "EBC-1",
    "CNhs11771": "stomach (fetal)",
    "CNhs12852": "submaxillary gland",
    "CNhs12318": "substantia nigra",
    "CNhs11244": "HS-SY-II",
    "CNhs11992": "synoviocytes",
    "CNhs12050": "synoviocytes",
    "CNhs10637": "temporal lobe",
    "CNhs11772": "temporal lobe (fetal)",
    "CNhs12997": "temporal lobe (fetal)",
    "CNhs12639": "tenocytes",
    "CNhs12640": "tenocytes",
    "CNhs12641": "tenocytes",
    "CNhs11878": "NCC-IT-A3",
    "CNhs11884": "NCR-G1",
    "CNhs11890": "PA-1",
    "CNhs11876": "ITO-II",
    "CNhs12351": "NEC14",
    "CNhs12362": "NEC15",
    "CNhs11726": "NEC8",
    "CNhs10632": "testis",
    "CNhs12998": "testis",
    "CNhs13794": "thalamus",
    "CNhs12314": "thalamus",
    "CNhs12858": "throat",
    "CNhs11770": "throat (fetal)",
    "CNhs10633": "thymus",
    "CNhs10650": "thymus (fetal)",
    "CNhs10634": "thyroid",
    "CNhs11872": "TCO-1",
    "CNhs11769": "thyroid (fetal)",
    "CNhs12853": "tongue",
    "CNhs11768": "tongue (fetal)",
    "CNhs10654": "tonsil",
    "CNhs11340": "trabecular meshwork cells",
    "CNhs12124": "trabecular meshwork cells",
    "CNhs10635": "trachea",
    "CNhs11766": "trachea (fetal)",
    "CNhs11092": "epithelial cells (trachea)",
    "CNhs11993": "epithelial cells (trachea)",
    "CNhs12051": "epithelial cells (trachea)",
    "CNhs10735": "5637",
    "CNhs11261": "JMSU1",
    "CNhs11828": "HGRT",
    "CNhs11883": "SUIT-2",
    "CNhs11765": "umbilical cord (fetal)",
    "CNhs10843": "urothelial cells",
    "CNhs11334": "urothelial cells",
    "CNhs12091": "urothelial cells",
    "CNhs12122": "urothelial cells",
    "CNhs11676": "uterus",
    "CNhs11763": "uterus (fetal)",
    "CNhs12854": "vagina",
    "CNhs12844": "vein",
    "CNhs11675": "blood (ribopure)",
    "CNhs11671": "blood (ribopure)",
    "CNhs11948": "blood (ribopure)",
    "CNhs11075": "blood (ribopure)",
    "CNhs11076": "blood (ribopure)",
    "CNhs11672": "blood (ribopure)",
    "CNhs11673": "blood (ribopure)",
    "CNhs11949": "blood (ribopure)",
    "CNhs11892": "G-401",
    "CNhs11728": "HFWT",
    "CNhs11813": "XPL 17",
}

#-------------#
# Classes     #
#-------------#

class Model(object): pass

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="this script inserts \"enhancer\" or \"tss\" data from FANTOM into GUD. download \"hg19_permissive_enhancers_expression_rle_tpm.csv.gz\" and \"hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz\" for enhancer and tss data, respectively.")

    parser.add_argument("matrix", help="Expression (TPM/RLE normalized) matrix across all FANTOM libraries")

    feats = ["enhancer", "tss"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature", metavar="feature_type")

    # Optional args
    parser.add_argument("-b", "--bed", help="BED file of features on which to focus (e.g. \"robust_enhancers.bed\")")
    parser.add_argument("-k", "--keep", action="store_true", help="Keep original sample fanmes from FANTOM (default = False)")
    parser.add_argument("--source", default="FANTOM", help="Source name (e.g. \"PMID:24670763\" for TSSs or \"PMID:24670764\" for enhancers; default = \"FANTOM\")")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default="hg19",
        help="Database name (default = \"hg19\")")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="Host name (default = localhost)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="Port number (default = 5506)")
    mysql_group.add_argument("-u", "--user", default=getpass.getuser(),
        help="User name (default = current user)")

    args = parser.parse_args()

    # Set default
    if not args.db:
        args.db = args.genome

    return args

def insert_fantom_to_gud_db(user, host, port, db, matrix_file,
    feat_type, source_name, bed_file=None, keep=False):

    # Initialize
    coordinates = set()
    fantom_sample_names = []
    db_name = "mysql://{}:@{}:{}/{}".format(
        user, host, port, db)
    if not database_exists(db_name):
        raise ValueError("GUD db does not exist: %s" % db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)
    today = str(date.today())
    if matrix_file.endswith(".gz"): gz = True
    else: gz = False

    # Initialize table
    if feat_type == "enhancer":
        table = Enhancer()
    if feat_type == "tss":
        table = TSS()
    if not engine.has_table(table.__tablename__):
        try:
            table.metadata.bind = engine
            table.metadata.create_all(engine)
        except:
            raise ValueError("Cannot create table: %s" % table.__tablename__)
    mapper(Model, table.__table__)

    # If BED file...
    if bed_file:
        # If BED file exists...
        if os.path.exists(bed_file):
            try:
                # Create BED object
                bed_obj = pybedtools.BedTool(bed_file)
                # Sort BED object
                for chrom, start, end in bed_obj.sort():
                    coordinates.add((chrom, int(start), int(end)))
            except:
                warnings.warn("\nCould not read file: \"%s\"\n\tSkipping file...\n" % file_name)

    # For each line...
    if feat_type == "enhancer":
        lines = GUDglobals.parse_csv_file(matrix_file, gz)
        counts_start_at = 1
    if feat_type == "tss":
        lines = GUDglobals.parse_tsv_file(matrix_file, gz)
        counts_start_at = 7
    for line in lines:
        # Skip comments
        if line[0].startswith("#"): continue
        # If no samples...
        if len(fantom_sample_names) == 0:
            for sample in line[counts_start_at:]:
                fantom_sample_names.append(unquote(sample))
        # ... Else...
        elif line[0].startswith("chr") or line[0].startswith("\"chr"):
            # Initialize
            samples = {}
            total_cages = 0.0
            # Get chrom, start, end
            if feat_type == "enhancer":
                m = re.search("(chr\S+)\:(\d+)\-(\d+)", line[0])
            if feat_type == "tss":
                m = re.search("(chr\S+)\:(\d+)\.\.(\d+),(\S)", line[0])
                n = re.search("p(\d+)@(\w+)", line[1])
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            if feat_type == "tss":
                strand = m.group(4)
            # Skip coordiantes
            if (chrom, start, end) in coordinates: continue
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\w{1,2})$", chrom)
            if not m.group(1) in GUDglobals.chroms: continue
            # Create model
            model = Model()
            model.bin = assign_bin(int(start), int(end))
            model.chrom = chrom
            model.start = start
            model.end = end
            model.experiment_type = "CAGE"
            model.source_name = source_name
            model.date = today
            if feat_type == "tss":
                model.gene = "n/a"
                model.tss = 1
                model.strand = strand
                if n:
                    model.gene = n.group(2)
                    model.tss = n.group(1)
            # For each sample...
            for i in range(counts_start_at, len(line)):
                # Initialize
                cages = "%.3f" % float(line[i])
                sample = fantom_sample_names[i - counts_start_at]
                # Keep original sample names
                if keep:
                    samples.setdefault(sample, [float(cages)])
                    total_cages += float(cages)
                else:
                    m = re.search("(CNhs\d+)", sample)
                    if m.group(1) in sample_names:
                        samples.setdefault(sample_names[m.group(1)], [])
                        samples[sample_names[m.group(1)]].append(float(cages))
                        total_cages += float(cages)
            # For each sample...
            for sample in samples:
                model.cell_or_tissue = sample
                if feat_type == "enhancer":
                    # Skip enhancers with 0 cages
                    if sum(samples[sample]) == 0:
                        continue
                    # Upsert model & commit
                    session.merge(model)
                    session.commit()
                if feat_type == "tss":
                    # For each id...
                    for i in range(len(samples[sample])):
                        model.replicate = i + 1
                        model.tpm = samples[sample][i]
                        model.percent_tpm = samples[sample][i] / total_cages
                        # Upsert model & commit
                        session.merge(model)
                        session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert FANTOM data to GUD database
    insert_fantom_to_gud_db(args.user, args.host, args.port,
        args.db, args.matrix, args.feat_type, args.source,
        args.bed, args.keep)