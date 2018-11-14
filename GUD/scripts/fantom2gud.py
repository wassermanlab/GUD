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

grouped_sample_names = {
    "acantholytic squamous carcinoma cell line:HCC1806 : CNhs1184": "HCC1806",
    "acute lymphoblastic leukemia (B-ALL) cell line:BALL-1 : CNhs1125": "BALL-1",
    "acute lymphoblastic leukemia (B-ALL) cell line:NALM-6 : CNhs1128": "NALM-6",
    "acute lymphoblastic leukemia (T-ALL) cell line:HPB-ALL : CNhs1074": "HPB-ALL",
    "acute lymphoblastic leukemia (T-ALL) cell line:Jurkat : CNhs1125": "Jurkat",
    "acute myeloid leukemia (FAB M0) cell line:KG-1 : CNhs1305": "KG-1",
    "acute myeloid leukemia (FAB M1) cell line:HYT-1 : CNhs1305": "HYT-1",
    "acute myeloid leukemia (FAB M2) cell line:Kasumi-1 : CNhs1350": "Kasumi-1",
    "acute myeloid leukemia (FAB M2) cell line:Kasumi-6 : CNhs1305": "Kasumi-6",
    "acute myeloid leukemia (FAB M2) cell line:NKM-1 : CNhs1186": "NKM-1",
    "acute myeloid leukemia (FAB M3) cell line:HL60 : CNhs1305": "HL60",
    "acute myeloid leukemia (FAB M4) cell line:FKH-1 : CNhs1350": "FKH-1",
    "acute myeloid leukemia (FAB M4) cell line:HNT-34 : CNhs1350": "HNT-34",
    "acute myeloid leukemia (FAB M4eo) cell line:EoL-1 : CNhs1305": "EoL-1",
    "acute myeloid leukemia (FAB M4eo) cell line:EoL-3 : CNhs1305": "EoL-3",
    "acute myeloid leukemia (FAB M5) cell line:NOMO-1 : CNhs1305": "NOMO-1",
    "acute myeloid leukemia (FAB M5) cell line:P31/FUJ : CNhs1305": "P31/FUJ",
    "acute myeloid leukemia (FAB M5) cell line:THP-1 (fresh) : CNhs1072": "THP-1",
#    "acute myeloid leukemia (FAB M5) cell line:THP-1 (revived) : CNhs1072": "THP-1",
#    "acute myeloid leukemia (FAB M5) cell line:THP-1 (thawed) : CNhs1072": "THP-1",
    "acute myeloid leukemia (FAB M5) cell line:U-937 DE-4 : CNhs1305": "U-937 DE-4",
    "acute myeloid leukemia (FAB M6) cell line:EEB : CNhs1305": "EEB",
    "acute myeloid leukemia (FAB M6) cell line:F-36E : CNhs1306": "F-36E",
    "acute myeloid leukemia (FAB M6) cell line:F-36P : CNhs1350": "F-36P",
    "acute myeloid leukemia (FAB M7) cell line:MKPL-1 : CNhs1188": "MKPL-1",
    "acute myeloid leukemia (FAB M7) cell line:M-MOK : CNhs1304": "M-MOK",
    "adenocarcinoma cell line:IM95m : CNhs1188": "IM95m",
    "Adipocyte - breast donor1 : CNhs1105": "adipocyte (breast)",
    "Adipocyte - breast donor2 : CNhs1196": "adipocyte (breast)",
    "Adipocyte - omental donor1 : CNhs1105": "adipocyte (omentum)",
    "Adipocyte - omental donor2 : CNhs1206": "adipocyte (omentum)",
    "Adipocyte - omental donor3 : CNhs1206": "adipocyte (omentum)",
    "Adipocyte - perirenal donor1 : CNhs1206": "adipocyte (perirenal)",
    "Adipocyte - subcutaneous donor1 : CNhs1249": "adipocyte (subcutaneous)",
    "Adipocyte - subcutaneous donor2 : CNhs1137": "adipocyte (subcutaneous)",
    "Adipocyte - subcutaneous donor3 : CNhs1201": "adipocyte (subcutaneous)",
    "adipose tissue adult pool1 : CNhs1061": "adipose tissue",
    "adrenal cortex adenocarcinoma cell line:SW-13 : CNhs1189": "SW-13",
    "adult T-cell leukemia cell line:ATN-1 : CNhs1073": "ATN-1",
    "alveolar cell carcinoma cell line:SW 1573 : CNhs1183": "SW 1573",
    "Alveolar Epithelial Cells donor1 : CNhs1132": "epithelial cells (alveolus)",
    "Alveolar Epithelial Cells donor2 : CNhs1208": "epithelial cells (alveolus)",
    "Amniotic Epithelial Cells donor1 : CNhs1134": "epithelial cells (amniotic membrane)",
    "Amniotic Epithelial Cells donor3 : CNhs1212": "epithelial cells (amniotic membrane)",
    "amniotic membrane cells donor1 : CNhs1250": "amniotic membrane cells",
    "amniotic membrane cells donor2 : CNhs1250": "amniotic membrane cells",
    "amniotic membrane cells donor3 : CNhs1237": "amniotic membrane cells",
    "amygdala - adult donor10196 : CNhs1379": "amygdala",
    "amygdala adult donor10252 : CNhs1231": "amygdala",
    "anaplastic carcinoma cell line:8305C : CNhs1074": "8305C",
    "anaplastic large cell lymphoma cell line:Ki-JK : CNhs1188": "Ki-JK",
    "anaplastic squamous cell carcinoma cell line:RPMI 2650 : CNhs1188": "RPMI 2650",
    "Anulus Pulposus Cell donor1 : CNhs1087": "anulus pulposus cells",
    "Anulus Pulposus Cell donor2 : CNhs1206": "anulus pulposus cells",
    "aorta adult pool1 : CNhs1176": "aorta",
    "appendix adult : CNhs1284": "appendix",
    "argyrophil small cell carcinoma cell line:TC-YIK : CNhs1172": "TC-YIK",
    "Astrocyte - cerebellum donor1 : CNhs1132": "astrocytes (cerebellum)",
    "Astrocyte - cerebellum donor2 : CNhs1208": "astrocytes (cerebellum)",
    "Astrocyte - cerebellum donor3 : CNhs1211": "astrocytes (cerebellum)",
    "Astrocyte - cerebral cortex donor1 : CNhs1086": "astrocytes (cerebral cortex)",
    "Astrocyte - cerebral cortex donor2 : CNhs1196": "astrocytes (cerebral cortex)",
    "Astrocyte - cerebral cortex donor3 : CNhs1200": "astrocytes (cerebral cortex)",
    "astrocytoma cell line:TM-31 : CNhs1074": "TM-31",
    "basal cell carcinoma cell line:TE 354.T : CNhs1193": "TE 354.T",
    "Basophils donor3 : CNhs1257": "basophils",
    "b cell line:RPMI1788 : CNhs1074": "RPMI1788",
    "bile duct carcinoma cell line:HuCCT1 : CNhs1075": "HuCCT1",
    "bile duct carcinoma cell line:TFK-1 : CNhs1126": "TFK-1",
    "biphenotypic B myelomonocytic leukemia cell line:MV-4-11 : CNhs1184": "MV-4-11",
    "bladder adult pool1 : CNhs1061": "bladder",
    "blood adult pool1 : CNhs1176": "blood",
    "B lymphoblastoid cell line: GM12878 ENCODE biol_rep1 : CNhs1233": "GM12878",
    "B lymphoblastoid cell line: GM12878 ENCODE biol_rep2 : CNhs1233": "GM12878",
    "B lymphoblastoid cell line: GM12878 ENCODE biol_rep3 : CNhs1233": "GM12878",
    "bone marrow stromal cell line:StromaNKtert : CNhs1193": "StromaNKtert",
    "brain adult donor1 : CNhs1179": "brain",
    "brain adult pool1 : CNhs1061": "brain",
    "brain fetal pool1 : CNhs1179": "brain (fetal)",
    "breast adult donor1 : CNhs1179": "breast",
    "breast carcinoma cell line:MCF7 : CNhs1194": "MCF-7",
    "breast carcinoma cell line:MDA-MB-453 : CNhs1073": "MDA-MB-453",
    "Bronchial Epithelial Cell donor4 : CNhs1205": "epithelial cells (bronchi)",
    "Bronchial Epithelial Cell donor5 : CNhs1205": "epithelial cells (bronchi)",
    "Bronchial Epithelial Cell donor6 : CNhs1206": "epithelial cells (bronchi)",
    "bronchial squamous cell carcinoma cell line:KNS-62 : CNhs1186": "KNS-62",
    "bronchioalveolar carcinoma cell line:NCI-H358 : CNhs1184": "NCI-H358",
    "bronchogenic carcinoma cell line:ChaGo-K-1 : CNhs1184": "ChaGo-K-1",
    "Burkitt's lymphoma cell line:DAUDI : CNhs1073": "DAUDI",
    "Burkitt's lymphoma cell line:RAJI : CNhs1126": "RAJI",
    "carcinoid cell line:NCI-H1770 : CNhs1183": "NCI-H1770",
    "carcinoid cell line:SK-PN-DW : CNhs1184": "SK-PN-DW",
    "carcinosarcoma cell line:JHUCS-1 : CNhs1174": "JHUCS-1",
    "Cardiac Myocyte donor1 : CNhs1234": "myocytes (heart)",
    "Cardiac Myocyte donor2 : CNhs1235": "myocytes (heart)",
    "Cardiac Myocyte donor3 : CNhs1257": "myocytes (heart)",
    "caudate nucleus - adult donor10196 : CNhs1380": "caudate nucleus",
    "caudate nucleus adult donor10252 : CNhs1232": "caudate nucleus",
    "CD14+CD16- Monocytes donor1 : CNhs1322": "monocytes (CD14-positive CD16-negative)",
    "CD14+CD16+ Monocytes donor1 : CNhs1354": "monocytes (CD14-positive CD16-positive)",
    "CD14-CD16+ Monocytes donor2 : CNhs1320": "monocytes (CD14-negative CD16-positive)",
    "CD14+CD16+ Monocytes donor2 : CNhs1320": "monocytes (CD14-positive CD16-positive)",
    "CD14+CD16- Monocytes donor2 : CNhs1321": "monocytes (CD14-positive CD16-negative)",
    "CD14-CD16+ Monocytes donor3 : CNhs1354": "monocytes (CD14-negative CD16-positive)",
    "CD14+CD16- Monocytes donor3 : CNhs1354": "monocytes (CD14-positive CD16-negative)",
    "CD14+CD16+ Monocytes donor3 : CNhs1354": "monocytes (CD14-positive CD16-positive)",
    "CD14+ monocyte derived endothelial progenitor cells donor1 : CNhs1085": "endothelial progenitor cells (derived from CD14-positive monocytes)",
    "CD14+ monocyte derived endothelial progenitor cells donor2 : CNhs1189": "endothelial progenitor cells (derived from CD14-positive monocytes)",
    "CD14+ monocyte derived endothelial progenitor cells donor3 : CNhs1190": "endothelial progenitor cells (derived from CD14-positive monocytes)",
    "CD14+ Monocytes donor1 : CNhs1085": "monocytes (CD14-positive)",
    "CD14+ Monocytes donor2 : CNhs1195": "monocytes (CD14-positive)",
    "CD14+ Monocytes donor3 : CNhs1199": "monocytes (CD14-positive)",
#    "CD14+ monocytes - mock treated donor1 : CNhs1346": "CD14-positive monocytes",
#    "CD14+ monocytes - mock treated donor2 : CNhs1348": "CD14-positive monocytes",
#    "CD14+ monocytes - mock treated donor3 : CNhs1349": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with BCG donor1 : CNhs1346": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with BCG donor2 : CNhs1347": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with BCG donor3 : CNhs1354": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with B-glucan donor1 : CNhs1347": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with B-glucan donor2 : CNhs1348": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with B-glucan donor3 : CNhs1349": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Candida donor1 : CNhs1347": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Candida donor2 : CNhs1348": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Candida donor3 : CNhs1349": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Cryptococcus donor1 : CNhs1347": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Cryptococcus donor2 : CNhs1348": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Cryptococcus donor3 : CNhs1354": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Group A streptococci donor1 : CNhs1346": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Group A streptococci donor2 : CNhs1353": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Group A streptococci donor3 : CNhs1349": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with IFN + N-hexane donor1 : CNhs1346": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with IFN + N-hexane donor2 : CNhs1347": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with IFN + N-hexane donor3 : CNhs1349": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with lipopolysaccharide donor1 : CNhs1347": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with lipopolysaccharide donor2 : CNhs1353": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with lipopolysaccharide donor3 : CNhs1354": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Salmonella donor1 : CNhs1347": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Salmonella donor2 : CNhs1348": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Salmonella donor3 : CNhs1349": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor1 : CNhs1346": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor2 : CNhs1348": "CD14-positive monocytes",
#    "CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor3 : CNhs1354": "CD14-positive monocytes",
    "CD19+ B Cells donor1 : CNhs1234": "B cells (CD19-positive)",
    "CD19+ B Cells donor2 : CNhs1235": "B cells (CD19-positive)",
    "CD19+ B Cells donor3 : CNhs1235": "B cells (CD19-positive)",
    "CD34+ stem cells - adult bone marrow derived donor1 tech_rep1 : CNhs1258": "CD34-positive stem cells (derived from bone marrow)",
    "CD4+CD25-CD45RA- memory conventional T cells donor3 : CNhs1353": "T cells (CD4-positive CD25-negative CD45RA-negative memory conventional)",
#    "CD4+CD25-CD45RA- memory conventional T cells expanded donor1 : CNhs1321": "CD4-positive CD25-negative CD45RA-negative memory conventional T cells",
    "CD4+CD25+CD45RA- memory regulatory T cells donor1 : CNhs1319": "T cells (CD4-positive CD25-positive CD45RA-negative memory regulatory)",
    "CD4+CD25+CD45RA- memory regulatory T cells donor2 : CNhs1320": "T cells (CD4-positive CD25-positive CD45RA-negative memory regulatory)",
    "CD4+CD25+CD45RA- memory regulatory T cells donor3 : CNhs1353": "T cells (CD4-positive CD25-positive CD45RA-negative memory regulatory)",
#    "CD4+CD25+CD45RA- memory regulatory T cells expanded donor1 : CNhs1320": "CD4-positive CD25-negative CD45RA-negative memory regulatory T cells",
#    "CD4+CD25+CD45RA- memory regulatory T cells expanded donor2 : CNhs1381": "CD4-positive CD25-negative CD45RA-negative memory regulatory T cells",
#    "CD4+CD25+CD45RA- memory regulatory T cells expanded donor3 : CNhs1381": "CD4-positive CD25-negative CD45RA-negative memory regulatory T cells",
    "CD4+CD25-CD45RA+ naive conventional T cells donor1 : CNhs1322": "T cells (CD4-positive CD25-negative CD45RA-positive naive conventional)",
    "CD4+CD25-CD45RA+ naive conventional T cells donor2 : CNhs1320": "T cells (CD4-positive CD25-negative CD45RA-positive naive conventional)",
    "CD4+CD25-CD45RA+ naive conventional T cells donor3 : CNhs1351": "T cells (CD4-positive CD25-negative CD45RA-positive naive conventional)",
#    "CD4+CD25-CD45RA+ naive conventional T cells expanded donor1 : CNhs1320": "CD4-positive CD25-negative CD45RA-positive naive conventional T cells",
#    "CD4+CD25-CD45RA+ naive conventional T cells expanded donor2 : CNhs1381": "CD4-positive CD25-negative CD45RA-positive naive conventional T cells",
#    "CD4+CD25-CD45RA+ naive conventional T cells expanded donor3 : CNhs1381": "CD4-positive CD25-negative CD45RA-positive naive conventional T cells",
    "CD4+CD25+CD45RA+ naive regulatory T cells donor3 : CNhs1351": "T cells (CD4-positive CD25-positive CD45RA-positive naive regulatory)",
#    "CD4+CD25+CD45RA+ naive regulatory T cells expanded donor1 : CNhs1320": "CD4-positive CD25-positive CD45RA-positive naive regulatory T cells",
    "CD4+ T Cells donor1 : CNhs1085": "T cells (CD4-positive)",
    "CD4+ T Cells donor2 : CNhs1195": "T cells (CD4-positive)",
    "CD4+ T Cells donor3 : CNhs1199": "T cells (CD4-positive)",
    "CD8+ T Cells donor1 : CNhs1085": "T cells (CD8-positive)",
    "CD8+ T Cells donor2 : CNhs1195": "T cells (CD8-positive)",
    "CD8+ T Cells donor3 : CNhs1199": "T cells (CD8-positive)",
    "cerebellum - adult donor10196 : CNhs1379": "cerebellum",
    "cerebellum adult donor10252 : CNhs1232": "cerebellum",
    "cerebellum adult pool1 : CNhs1179": "cerebellum",
    "cerebral meninges adult : CNhs1284": "meninges (brain)",
    "cervical cancer cell line:D98-AH2 : CNhs1128": "D98-AH2",
    "cervical cancer cell line:ME-180 : CNhs1128": "ME-180",
    "cervix adult pool1 : CNhs1061": "cervix",
    "cholangiocellular carcinoma cell line:HuH-28 : CNhs1128": "HuH-28",
#    "Chondrocyte - de diff donor1 : CNhs1192": ,
#    "Chondrocyte - de diff donor2 : CNhs1137": ,
#    "Chondrocyte - de diff donor3 : CNhs1202": ,
#    "Chondrocyte - re diff donor2 : CNhs1137": ,
#    "Chondrocyte - re diff donor3 : CNhs1202": ,
    "choriocarcinoma cell line:BeWo : CNhs1074": "BeWo",
    "choriocarcinoma cell line:SCH : CNhs1187": "SCH",
    "choriocarcinoma  cell line:T3M-3 : CNhs1182": "T3M-3",
    "chorionic membrane cells donor1 : CNhs1250": "chorionic membrane cells",
    "chorionic membrane cells donor2 : CNhs1250": "chorionic membrane cells",
    "chorionic membrane cells donor3 : CNhs1238": "chorionic membrane cells",
    "chronic lymphocytic leukemia (T-CLL) cell line:SKW-3 : CNhs1171": "SKW-3",
    "chronic myeloblastic leukemia (CML) cell line:KCL-22 : CNhs1188": "KCL-22",
    "chronic myelogenous leukemia cell line:K562 : CNhs1125": "K562",
    "chronic myelogenous leukemia cell line:K562 ENCODE biol_rep1 : CNhs1233": "K562",
    "chronic myelogenous leukemia cell line:K562 ENCODE biol_rep2 : CNhs1233": "K562",
    "chronic myelogenous leukemia cell line:K562 ENCODE biol_rep3 : CNhs1233": "K562",
    "chronic myelogenous leukemia cell line:KU812 : CNhs1072": "KU812",
    "chronic myelogenous leukemia (CML) cell line:MEG-A2 : CNhs1186": "MEG-A2",
    "Ciliary Epithelial Cells donor1 : CNhs1087": "ciliary epithelial cells",
    "Ciliary Epithelial Cells donor2 : CNhs1196": "ciliary epithelial cells",
    "Ciliary Epithelial Cells donor3 : CNhs1200": "ciliary epithelial cells",
    "clear cell carcinoma cell line:JHOC-5 : CNhs1174": "JHOC-5",
    "clear cell carcinoma cell line:TEN : CNhs1193": "TEN",
#    "Clontech Human Universal Reference Total RNA pool1 : CNhs1060": ,
    "colon adult donor1 : CNhs1179": "colon",
    "colon adult pool1 : CNhs1061": "colon",
    "colon carcinoma cell line:CACO-2 : CNhs1128": "CACO-2",
    "colon carcinoma cell line:COLO-320 : CNhs1073": "COLO-320",
    "colon fetal donor1 : CNhs1178": "colon (fetal)",
#    "cord blood derived cell line:COBL-a 24h infection(-C) : CNhs1104": ,
#    "cord blood derived cell line:COBL-a 24h infection : CNhs1105": ,
    "cord blood derived cell line:COBL-a untreated : CNhs1104": "COBL-a",
    "Corneal Epithelial Cells donor1 : CNhs1133": "epithelial cells (cornea)",
    "Corneal Epithelial Cells donor3 : CNhs1212": "epithelial cells (cornea)",
    "corpus callosum adult pool1 : CNhs1064": "corpus callosum",
    "Dendritic Cells - monocyte immature derived donor1 tech_rep1 : CNhs1085": "dendritic cells (derived from immature monocytes)",
    "Dendritic Cells - monocyte immature derived donor1 tech_rep2 : CNhs1106": "dendritic cells (derived from immature monocytes)",
    "Dendritic Cells - monocyte immature derived donor3 : CNhs1200": "dendritic cells (derived from immature monocytes)",
    "Dendritic Cells - plasmacytoid donor1 : CNhs1085": "plasmacytoid dendritic cells",
    "diaphragm fetal donor1 : CNhs1177": "diaphragm (fetal)",
    "diencephalon adult : CNhs1261": "diencephalon",
    "diffuse large B-cell lymphoma cell line:CTB-1 : CNhs1174": "CTB-1",
    "ductal cell carcinoma cell line:KLM-1 : CNhs1110": "KLM-1",
    "ductal cell carcinoma cell line:MIA Paca2 : CNhs1125": "MIA Paca2",
    "ductus deferens adult : CNhs1284": "ductus deferens",
    "duodenum fetal donor1 tech_rep1 : CNhs1178": "duodenum (fetal)",
    "duodenum fetal donor1 tech_rep2 : CNhs1299": "duodenum (fetal)",
    "dura mater adult donor1 : CNhs1064": "dura mater",
#    "embryonic kidney cell line: HEK293/SLAM infection 24hr : CNhs1104": "HEK293/SLAM",
    "embryonic kidney cell line: HEK293/SLAM untreated : CNhs1104": "HEK293/SLAM",
    "embryonic pancreas cell line:1B2C6 : CNhs1173": "1B2C6",
    "embryonic pancreas cell line:1C3D3 : CNhs1173": "1C3D3",
    "embryonic pancreas cell line:1C3IKEI : CNhs1173": "1C3IKEI",
    "embryonic pancreas cell line:2C6 : CNhs1181": "2C6",
    "endometrial carcinoma cell line:OMC-2 : CNhs1126": "OMC-2",
    "endometrial stromal sarcoma cell line:OMC-9 : CNhs1124": "OMC-9",
    "endometrioid adenocarcinoma cell line:JHUEM-1 : CNhs1174": "JHUEM-1",
    "Endothelial Cells - Aortic donor0 : CNhs1083": "endothelial cells (aorta)",
    "Endothelial Cells - Aortic donor1 : CNhs1249": "endothelial cells (aorta)",
    "Endothelial Cells - Aortic donor2 : CNhs1137": "endothelial cells (aorta)",
    "Endothelial Cells - Aortic donor3 : CNhs1202": "endothelial cells (aorta)",
    "Endothelial Cells - Artery donor1 : CNhs1249": "endothelial cells (artery)",
    "Endothelial Cells - Artery donor2 : CNhs1197": "endothelial cells (artery)",
    "Endothelial Cells - Artery donor3 : CNhs1202": "endothelial cells (artery)",
    "Endothelial Cells - Lymphatic donor1 : CNhs1086": "endothelial cells (lymph node)",
    "Endothelial Cells - Lymphatic donor2 : CNhs1190": "endothelial cells (lymph node)",
    "Endothelial Cells - Lymphatic donor3 : CNhs1190": "endothelial cells (lymph )",
    "Endothelial Cells - Microvascular donor1 : CNhs1192": "endothelial cells (microvasculature)",
    "Endothelial Cells - Microvascular donor2 : CNhs1137": "endothelial cells (microvasculature)",
    "Endothelial Cells - Microvascular donor3 : CNhs1202": "endothelial cells (microvasculature)",
    "Endothelial Cells - Thoracic donor1 : CNhs1192": "endothelial cells (thorax)",
    "Endothelial Cells - Thoracic donor2 : CNhs1197": "endothelial cells (thorax)",
    "Endothelial Cells - Umbilical vein donor1 : CNhs1087": "endothelial cells (umbilical vein)",
    "Endothelial Cells - Umbilical vein donor2 : CNhs1196": "endothelial cells (umbilical vein)",
    "Endothelial Cells - Umbilical vein donor3 : CNhs1201": "endothelial cells (umbilical vein)",
    "Endothelial Cells - Vein donor1 : CNhs1249": "endothelial cells (vein)",
    "Endothelial Cells - Vein donor2 : CNhs1137": "endothelial cells (vein)",
    "Endothelial Cells - Vein donor3 : CNhs1202": "endothelial cells (vein)",
    "epidermoid carcinoma  cell line:A431 : CNhs1074": "A431",
    "epidermoid carcinoma cell line:Ca Ski : CNhs1074": "Ca Ski",
    "epididymis adult : CNhs1284": "epididymis",
    "epithelioid sarcoma cell line:HS-ES-1 : CNhs1124": "HS-ES-1",
    "epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep1 : CNhs1232": "HeLa-S3",
    "epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep2 : CNhs1232": "HeLa-S3",
    "epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep3 : CNhs1232": "HeLa-S3",
    "Esophageal Epithelial Cells donor1 : CNhs1132": "epithelial cells (esophagus)",
    "esophagus adult pool1 : CNhs1062": "esophagus",
    "Ewing's sarcoma cell line:Hs 863.T : CNhs1183": "Hs 863.T",
    "extraskeletal myxoid chondrosarcoma cell line:H-EMC-SS : CNhs1072": "H-EMC-SS",
    "eye fetal donor1 : CNhs1176": "eye (fetal)",
    "Fibroblast - Aortic Adventitial donor1 : CNhs1087": "fibroblast (aortic adventitia)",
    "Fibroblast - Aortic Adventitial donor3 : CNhs1201": "fibroblast (aortic adventitia)",
    "Fibroblast - Cardiac donor1 : CNhs1249": "fibroblast (heart)",
    "Fibroblast - Cardiac donor2 : CNhs1137": "fibroblast (heart)",
    "Fibroblast - Cardiac donor3 : CNhs1202": "fibroblast (heart)",
    "Fibroblast - Cardiac donor4 : CNhs1190": "fibroblast (heart)",
    "Fibroblast - Cardiac donor5 : CNhs1205": "fibroblast (heart)",
    "Fibroblast - Cardiac donor6 : CNhs1206": "fibroblast (heart)",
    "Fibroblast - Choroid Plexus donor1 : CNhs1131": "fibroblast (choroid plexus)",
    "Fibroblast - Choroid Plexus donor2 : CNhs1234": "fibroblast (choroid plexus)",
    "Fibroblast - Conjunctival donor1 : CNhs1133": "fibroblast (conjunctiva)",
    "Fibroblast - Dermal donor1 : CNhs1249": "fibroblast (dermis)",
    "Fibroblast - Dermal donor2 : CNhs1137": "fibroblast (dermis)",
    "Fibroblast - Dermal donor3 : CNhs1202": "fibroblast (dermis)",
    "Fibroblast - Dermal donor4 : CNhs1205": "fibroblast (dermis)",
    "Fibroblast - Dermal donor5 : CNhs1205": "fibroblast (dermis)",
    "Fibroblast - Dermal donor6 : CNhs1205": "fibroblast (dermis)",
    "Fibroblast - Gingival donor1 : CNhs1086": "fibroblast (gingiva)",
    "Fibroblast - Gingival donor2 : CNhs1196": "fibroblast (gingiva)",
    "Fibroblast - Gingival donor3 : CNhs1200": "fibroblast (gingiva)",
    "Fibroblast - Gingival donor4 (GFH2) : CNhs1084": "fibroblast (gingiva)",
    "Fibroblast - Gingival donor5 (GFH3) : CNhs1195": "fibroblast (gingiva)",
    "Fibroblast - Lymphatic donor1 : CNhs1132": "fibroblast (lymph node)",
    "Fibroblast - Lymphatic donor3 : CNhs1211": "fibroblast (lymph node)",
    "Fibroblast - Periodontal Ligament donor1 : CNhs1086": "fibroblast (periodontal ligament)",
    "Fibroblast - Periodontal Ligament donor2 : CNhs1196": "fibroblast (periodontal ligament)",
    "Fibroblast - Periodontal Ligament donor3 : CNhs1190": "fibroblast (periodontal ligament)",
    "Fibroblast - Periodontal Ligament donor4 (PL29) : CNhs1249": "fibroblast (periodontal ligament)",
    "Fibroblast - Periodontal Ligament donor5 (PL30) : CNhs1195": "fibroblast (periodontal ligament)",
    "Fibroblast - Periodontal Ligament donor6 (PLH3) : CNhs1199": "fibroblast (periodontal ligament)",
    "Fibroblast - Pulmonary Artery donor1 : CNhs1087": "fibroblast (pulmonary artery)",
#    "Fibroblast - skin dystrophia myotonica donor1 : CNhs1135": "fibroblast (skin of patient with dystrophia myotonica)",
#    "Fibroblast - skin dystrophia myotonica donor2 : CNhs1135": "fibroblast (skin of patient with dystrophia myotonica)",
#    "Fibroblast - skin dystrophia myotonica donor3 : CNhs1191": "fibroblast (skin of patient with dystrophia myotonica)",
    "Fibroblast - skin normal donor1 : CNhs1135": "fibroblast (skin)",
    "Fibroblast - skin normal donor2 : CNhs1191": "fibroblast (skin)",
#    "Fibroblast - skin spinal muscular atrophy donor1 : CNhs1107": "fibroblast (skin of patient with spinal muscular atrophy)",
#    "Fibroblast - skin spinal muscular atrophy donor2 : CNhs1191": "fibroblast (skin of patient with spinal muscular atrophy)",
#    "Fibroblast - skin spinal muscular atrophy donor3 : CNhs1191": "fibroblast (skin of patient with spinal muscular atrophy)",
#    "Fibroblast - skin walker warburg donor1 : CNhs1135": "fibroblast (skin of patient with Walker-Warburg syndrome)",
    "fibrosarcoma cell line:HT-1080 : CNhs1186": "HT-1080",
    "fibrous histiocytoma cell line:GCT TIB-223 : CNhs1184": "GCT TIB-223",
    "frontal lobe adult pool1 : CNhs1064": "frontal lobe",
    "gall bladder adult : CNhs1284": "gall bladder",
    "gall bladder carcinoma cell line:TGBC14TKB : CNhs1125": "TGBC14TKB",
    "gall bladder carcinoma cell line:TGBC2TKB : CNhs1073": "TGBC2TKB",
    "gastric adenocarcinoma cell line:MKN1 : CNhs1173": "MKN1",
    "gastric adenocarcinoma cell line:MKN45 : CNhs1181": "MKN45",
    "gastric cancer cell line:AZ521 : CNhs1128": "AZ521",
    "gastrointestinal carcinoma cell line:ECC12 : CNhs1173": "ECC12",
    "giant cell carcinoma cell line:LU65 : CNhs1127": "LU65",
    "giant cell carcinoma cell line:Lu99B : CNhs1075": "Lu99B",
    "Gingival epithelial cells donor1 (GEA11) : CNhs1106": "epithelial cells (gingiva)",
    "Gingival epithelial cells donor2 (GEA14) : CNhs1189": "epithelial cells (gingiva)",
    "Gingival epithelial cells donor3 (GEA15) : CNhs1190": "epithelial cells (gingiva)",
    "glassy cell carcinoma cell line:HOKUG : CNhs1182": "HOKUG",
    "glioblastoma cell line:A172 : CNhs1118": "A172",
    "glioblastoma cell line:A172 tech_rep2 : CNhs1124": "A172",
    "glioblastoma cell line:T98G : CNhs1127": "T98G",
    "glioma cell line:GI-1 : CNhs1073": "GI-1",
    "globus pallidus - adult donor10196 : CNhs1380": "globus pallidus",
    "globus pallidus adult donor10252 : CNhs1231": "globus pallidus",
    "granulosa cell tumor cell line:KGN : CNhs1174": "KGN",
    "Hair Follicle Dermal Papilla Cells donor1 : CNhs1250": "hair follicle dermal papilla cells",
    "Hair Follicle Dermal Papilla Cells donor2 : CNhs1197": "hair follicle dermal papilla cells",
    "Hair Follicle Dermal Papilla Cells donor3 : CNhs1203": "hair follicle dermal papilla cells",
    "Hair Follicle Outer Root Sheath Cells donor1 : CNhs1233": "hair follicle outer root sheath cells",
    "Hair Follicle Outer Root Sheath Cells donor2 : CNhs1234": "hair follicle outer root sheath cells",
    "hairy cell leukemia cell line:Mo : CNhs1184": "Mo",
#    "heart adult diseased donor1 : CNhs1175": "heart from diseased patient",
#    "heart adult diseased post-infarction donor1 : CNhs1175": "heart from diseased patient (post-infarction)",
    "heart adult pool1 : CNhs1062": "heart",
    "heart fetal pool1 : CNhs1065": "heart (fetal)",
    "heart - mitral valve adult : CNhs1285": "mitral valve",
    "heart - pulmonic valve adult : CNhs1285": "pulmonic valve",
    "heart - tricuspid valve adult : CNhs1285": "tricuspid valve",
    "Hep-2 cells mock treated biol_rep1 : CNhs1347": "HEp-2",
    "Hep-2 cells mock treated biol_rep2 : CNhs1350": "HEp-2",
    "Hep-2 cells mock treated biol_rep3 : CNhs1350": "HEp-2",
#    "Hep-2 cells treated with Streptococci strain 5448 biol_rep1 : CNhs1347": "Hep-2",
#    "Hep-2 cells treated with Streptococci strain 5448 biol_rep2 : CNhs1349": "Hep-2",
#    "Hep-2 cells treated with Streptococci strain 5448 biol_rep3 : CNhs1349": "Hep-2",
#    "Hep-2 cells treated with Streptococci strain JRS4 biol_rep1 : CNhs1347": "Hep-2",
#    "Hep-2 cells treated with Streptococci strain JRS4 biol_rep2 : CNhs1349": "Hep-2",
#    "Hep-2 cells treated with Streptococci strain JRS4 biol_rep3 : CNhs1349": "Hep-2",
    "hepatic mesenchymal tumor cell line:LI90 : CNhs1186": "LI90",
    "Hepatic Sinusoidal Endothelial Cells donor1 : CNhs1207": "endothelial cells (liver sinusoid)",
    "Hepatic Sinusoidal Endothelial Cells donor2 : CNhs1209": "endothelial cells (liver sinusoid)",
    "Hepatic Stellate Cells (lipocyte) donor1 : CNhs1133": "hepatic stellate cells",
    "Hepatic Stellate Cells (lipocyte) donor2 : CNhs1209": "hepatic stellate cells",
    "hepatoblastoma cell line:HuH-6 : CNhs1174": "HuH-6",
    "hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep1 : CNhs1232": "HepG2",
    "hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep2 : CNhs1232": "HepG2",
    "hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep3 : CNhs1233": "HepG2",
    "Hepatocyte donor1 : CNhs1234": "hepatocytes",
    "Hepatocyte donor2 : CNhs1234": "hepatocytes",
    "Hepatocyte donor3 : CNhs1262": "hepatocytes",
    "hepatoma cell line:Li-7 : CNhs1127": "Li-7",
    "hereditary spherocytic anemia cell line:WIL2-NS : CNhs1189": "WIL2-NS",
    "hippocampus - adult donor10196 : CNhs1379": "hippocampus",
    "hippocampus adult donor10252 : CNhs1231": "hippocampus",
    "Hodgkin's lymphoma cell line:HD-Mar2 : CNhs1171": "HD-Mar2",
    "immature langerhans cells donor1 : CNhs1353": "langerhans cells (immature)",
    "immature langerhans cells donor2 : CNhs1348": "langerhans cells (immature)",
    "insula adult pool1 : CNhs1064": "insula",
    "Intestinal epithelial cells (polarized) donor1 : CNhs1087": "epithelial cells (intestine)",
    "Iris Pigment Epithelial Cells donor1 : CNhs1259": "iris pigment epithelium",
    "Keratinocyte - epidermal donor1 : CNhs1106": "keratinocytes (epidermis)",
    "Keratinocyte - epidermal donor2 : CNhs1138": "keratinocytes (epidermis)",
    "Keratinocyte - epidermal donor3 : CNhs1203": "keratinocytes (epidermis)",
    "Keratinocyte - oral donor1 : CNhs1087": "keratinocytes (oral)",
    "keratoacanthoma cell line:HKA-1 : CNhs1188": "HKA-1",
    "Keratocytes donor1 : CNhs1133": "keratocytes",
    "Keratocytes donor2 : CNhs1209": "keratocytes",
    "kidney adult pool1 : CNhs1062": "kidney",
    "kidney fetal pool1 : CNhs1065": "kidney (fetal)",
    "large cell lung carcinoma cell line:IA-LM : CNhs1127": "IA-LM",
    "large cell lung carcinoma cell line:NCI-H460 : CNhs1280": "NCI-H460",
    "large cell non-keratinizing squamous carcinoma cell line:SKG-II-SF : CNhs1182": "SKG-II-SF",
    "left atrium adult donor1 : CNhs1179": "left atrium",
    "left ventricle adult donor1 : CNhs1178": "left ventricle",
    "leiomyoblastoma  cell line:G-402 : CNhs1184": "G-402",
    "leiomyoma cell line:10964C : CNhs1172": "10964C",
    "leiomyoma cell line:15242A : CNhs1172": "15242A",
    "leiomyoma cell line:15425 : CNhs1172": "15425",
    "lens epithelial cell line:SRA 01/04 : CNhs1175": "SRA 01/04",
    "Lens Epithelial Cells donor1 : CNhs1234": "epithelial cells (lens)",
    "Lens Epithelial Cells donor2 : CNhs1256": "epithelial cells (lens)",
    "Lens Epithelial Cells donor3 : CNhs1257": "epithelial cells (lens)",
    "leukemia chronic megakaryoblastic  cell line:MEG-01 : CNhs1185": "MEG-01",
    "liposarcoma  cell line:KMLS-1 : CNhs1187": "KMLS-1",
    "liposarcoma  cell line:SW 872 : CNhs1185": "SW 872",
    "liver adult pool1 : CNhs1062": "liver",
    "liver fetal pool1 : CNhs1179": "liver (fetal)",
    "locus coeruleus - adult donor10196 : CNhs1380": "locus coeruleus",
    "locus coeruleus adult donor10252 : CNhs1232": "locus coeruleus",
    "lung adenocarcinoma cell line:A549 : CNhs1127": "A549",
    "lung adenocarcinoma cell line:PC-14 : CNhs1072": "PC-14",
    "lung adult pool1 : CNhs1062": "lung",
    "lung fetal donor1 : CNhs1168": "lung (fetal)",
    "lung right lower lobe adult donor1 : CNhs1178": "lung (right lower lobe)",
    "lymphangiectasia cell line:DS-1 : CNhs1185": "DS-1",
    "lymph node adult donor1 : CNhs1178": "lymph node",
    "lymphoma malignant hairy B-cell cell line:MLMA : CNhs1193": "MLMA",
    "Macrophage - monocyte derived donor1 : CNhs1086": "macrophages (derived from monocytes)",
    "Macrophage - monocyte derived donor2 : CNhs1189": "macrophages (derived from monocytes)",
    "Macrophage - monocyte derived donor3 : CNhs1200": "macrophages (derived from monocytes)",
    "malignant trichilemmal cyst cell line:DJM-1 : CNhs1073": "DJM-1",
#    "Mallassez-derived cells donor2 : CNhs1355": "mallassez epithelial cells",
#    "Mallassez-derived cells donor3 : CNhs1355": "mallassez epithelial cells",
    "Mammary Epithelial Cell donor1 : CNhs1107": "epithelial cells (breast)",
    "Mammary Epithelial Cell donor2 : CNhs1138": "epithelial cells (breast)",
    "Mammary Epithelial Cell donor3 : CNhs1203": "epithelial cells (breast)",
    "Mast cell donor1 : CNhs1256": "mast cells",
    "Mast cell donor2 : CNhs1259": "mast cells",
    "Mast cell donor3 : CNhs1259": "mast cells",
    "Mast cell donor4 : CNhs1259": "mast cells",
#    "Mast cell - stimulated donor1 : CNhs1107": "mast cells",
    "maxillary sinus tumor cell line:HSQ-89 : CNhs1073": "HSQ-89",
    "medial frontal gyrus - adult donor10196 : CNhs1379": "medial frontal gyrus",
    "medial temporal gyrus - adult donor10196 : CNhs1380": "medial temporal gyrus",
    "medial temporal gyrus adult donor10252 : CNhs1231": "medial temporal gyrus",
    "medulla oblongata - adult donor10196 : CNhs1380": "medulla oblongata",
    "medulla oblongata adult donor10252 : CNhs1231": "medulla oblongata",
    "medulla oblongata adult pool1 : CNhs1064": "medulla oblongata",
    "medulloblastoma  cell line:D283 Med : CNhs1280": "D283 Med",
    "medulloblastoma  cell line:ONS-76 : CNhs1186": "ONS-76",
    "Melanocyte - dark donor3 : CNhs1257": "melanocytes (dark donor)",
    "Melanocyte - light donor1 : CNhs1130": "melanocytes (light donor)",
    "Melanocyte - light donor2 : CNhs1138": "melanocytes (light donor)",
    "Melanocyte - light donor3 : CNhs1203": "melanocytes (light donor)",
    "melanoma cell line:COLO 679 : CNhs1128": "COLO 679",
    "melanoma cell line:G-361 : CNhs1125": "G-361",
    "Meningeal Cells donor1 : CNhs1132": "meningeal cells",
    "Meningeal Cells donor2 : CNhs1208": "meningeal cells",
    "Meningeal Cells donor3 : CNhs1273": "meningeal cells",
    "meningioma cell line:HKBMM : CNhs1194": "HKBMM",
    "merkel cell carcinoma cell line:MKL-1 : CNhs1283": "MKL-1",
    "merkel cell carcinoma cell line:MS-1 : CNhs1283": "MS-1",
    "mesenchymal precursor cell - adipose donor1 : CNhs1236": "mesenchymal precursor cells (adipose tissue)",
    "mesenchymal precursor cell - adipose donor2 : CNhs1236": "mesenchymal precursor cells (adipose tissue)",
    "mesenchymal precursor cell - adipose donor3 : CNhs1236": "mesenchymal precursor cells (adipose tissue)",
    "mesenchymal precursor cell - bone marrow donor1 : CNhs1236": "mesenchymal precursor cells (bone marrow)",
    "mesenchymal precursor cell - bone marrow donor2 : CNhs1236": "mesenchymal precursor cells (bone marrow)",
    "mesenchymal precursor cell - bone marrow donor3 : CNhs1309": "mesenchymal precursor cells (bone marrow)",
    "mesenchymal precursor cell - cardiac donor1 : CNhs1236": "mesenchymal precursor cells (heart)",
    "mesenchymal precursor cell - cardiac donor2 : CNhs1236": "mesenchymal precursor cells (heart)",
    "mesenchymal precursor cell - cardiac donor3 : CNhs1237": "mesenchymal precursor cells (heart)",
    "mesenchymal precursor cell - cardiac donor4 : CNhs1237": "mesenchymal precursor cells (heart)",
#    "mesenchymal precursor cell - ovarian cancer left ovary donor1 : CNhs1237": "mesenchymal precursor cells of left ovarian cancer",
#    "mesenchymal precursor cell - ovarian cancer left ovary donor2 : CNhs1309": "mesenchymal precursor cells of left ovarian cancer",
#    "mesenchymal precursor cell - ovarian cancer left ovary donor3 : CNhs1237": "mesenchymal precursor cells of left ovarian cancer",
#    "mesenchymal precursor cell - ovarian cancer left ovary donor4 : CNhs1309": "mesenchymal precursor cells of left ovarian cancer",
#    "mesenchymal precursor cell - ovarian cancer metastasis donor1 : CNhs1237": "mesenchymal precursor cells of ovarian cancer metastasis",
#    "mesenchymal precursor cell - ovarian cancer metastasis donor2 : CNhs1309": "mesenchymal precursor cells of ovarian cancer metastasis",
#    "mesenchymal precursor cell - ovarian cancer metastasis donor3 : CNhs1237": "mesenchymal precursor cells of ovarian cancer metastasis",
#    "mesenchymal precursor cell - ovarian cancer metastasis donor4 : CNhs1309": "mesenchymal precursor cells of ovarian cancer metastasis",
#    "mesenchymal precursor cell - ovarian cancer right ovary donor1 : CNhs1237": "mesenchymal precursor cells of right ovarian cancer",
#    "mesenchymal precursor cell - ovarian cancer right ovary donor2 : CNhs1237": "mesenchymal precursor cells of right ovarian cancer",
#    "mesenchymal precursor cell - ovarian cancer right ovary donor3 (SOC-57-02) : CNhs1237": "mesenchymal precursor cells of right ovarian cancer",
#    "mesenchymal precursor cell - ovarian cancer right ovary donor3 (SOC-57-02-G) : CNhs1350": "mesenchymal precursor cells of right ovarian cancer",
#    "mesenchymal precursor cell - ovarian cancer right ovary donor4 : CNhs1309": "mesenchymal precursor cells of right ovarian cancer",
    "mesenchymal stem cell line:Hu5/E18 : CNhs1171": "Hu5/E18",
    "Mesenchymal stem cells - adipose donor0 : CNhs1084": "mesenchymal stem cells (adipose tissue)",
    "Mesenchymal Stem Cells - adipose donor1 : CNhs1134": "mesenchymal stem cells (adipose tissue)",
    "Mesenchymal Stem Cells - adipose donor3 : CNhs1292": "mesenchymal stem cells (adipose tissue)",
    "Mesenchymal Stem Cells - amniotic membrane donor1 : CNhs1134": "mesenchymal stem cells (amniotic membrane)",
    "Mesenchymal Stem Cells - amniotic membrane donor2 : CNhs1210": "mesenchymal stem cells (amniotic membrane)",
    "Mesenchymal Stem Cells - bone marrow donor1 : CNhs1134": "mesenchymal stem cells (bone marrow)",
    "Mesenchymal Stem Cells - bone marrow donor2 : CNhs1210": "mesenchymal stem cells (bone marrow)",
    "Mesenchymal Stem Cells - bone marrow donor3 : CNhs1212": "mesenchymal stem cells (bone marrow)",
    "Mesenchymal stem cells - hepatic donor0 : CNhs1084": "mesenchymal stem cells (liver)",
    "Mesenchymal Stem Cells - hepatic donor2 : CNhs1273": "mesenchymal stem cells (liver)",
    "Mesenchymal stem cells - umbilical donor0 : CNhs1249": "mesenchymal stem cells (umbilical cord)",
    "Mesenchymal Stem Cells - umbilical donor1 : CNhs1134": "mesenchymal stem cells (umbilical cord)",
    "Mesenchymal Stem Cells - umbilical donor3 : CNhs1212": "mesenchymal stem cells (umbilical cord)",
#    "Mesenchymal Stem Cells - Wharton's Jelly donor1 : CNhs1105": "mesenchymal stem cells from Wharton's Jelly patient",
    "Mesothelial Cells donor1 : CNhs1085": "mesothelial cells",
    "Mesothelial Cells donor3 : CNhs1201": "mesothelial cells",
    "mesothelioma cell line:ACC-MESO-1 : CNhs1126": "ACC-MESO-1",
    "mesothelioma cell line:ACC-MESO-4 : CNhs1126": "ACC-MESO-4",
    "mesothelioma cell line:Mero-25 : CNhs1306": "Mero-25",
    "mesothelioma cell line:Mero-41 : CNhs1306": "Mero-41",
    "mesothelioma cell line:Mero-48a : CNhs1306": "Mero-48a",
    "mesothelioma cell line:Mero-82 : CNhs1306": "Mero-82",
    "mesothelioma cell line:Mero-83 : CNhs1307": "Mero-83",
    "mesothelioma cell line:Mero-84 : CNhs1307": "Mero-84",
    "mesothelioma cell line:Mero-95 : CNhs1307": "Mero-95",
    "mesothelioma cell line:NCI-H2052 : CNhs1306": "NCI-H2052",
    "mesothelioma cell line:NCI-H226 : CNhs1306": "NCI-H226",
    "mesothelioma cell line:NCI-H2452 : CNhs1306": "NCI-H2452",
    "mesothelioma cell line:NCI-H28 : CNhs1306": "NCI-H28",
    "mesothelioma cell line:No36 : CNhs1307": "No36",
    "mesothelioma cell line:ONE58 : CNhs1307": "ONE58",
    "middle temporal gyrus donor10252 : CNhs1231": "middle temporal gyrus",
    "migratory langerhans cells donor1 : CNhs1353": "langerhans cells (migratory)",
    "migratory langerhans cells donor2 : CNhs1353": "langerhans cells (migratory)",
    "migratory langerhans cells donor3 : CNhs1354": "langerhans cells (migratory)",
    "mixed mullerian tumor cell line:HTMMT : CNhs1194": "HTMMT",
    "mucinous adenocarcinoma cell line:JHOM-1 : CNhs1175": "JHOM-1",
    "mucinous cystadenocarcinoma  cell line:MCAS : CNhs1187": "MCAS",
    "Multipotent Cord Blood Unrestricted Somatic Stem Cells donor1 : CNhs1135": "multipotent stem cells (umbilical cord blood)",
    "Multipotent Cord Blood Unrestricted Somatic Stem Cells donor2 : CNhs1210": "multipotent stem cells (umbilical cord blood)",
    "mycosis fungoides T cell lymphoma cell line:HuT 102 TIB-162 : CNhs1185": "HuT 102 TIB-162",
    "myelodysplastic syndrome cell line:SKM-1 : CNhs1193": "SKM-1",
    "myeloma cell line:PCM6 : CNhs1125": "PCM6",
    "Myoblast donor1 : CNhs1087": "myoblasts",
    "Myoblast donor2 : CNhs1196": "myoblasts",
    "Myoblast donor3 : CNhs1190": "myoblasts",
    "myxofibrosarcoma cell line:MFH-ino : CNhs1172": "MFH-ino",
    "myxofibrosarcoma cell line:NMFH-1 : CNhs1182": "NMFH-1",
    "nasal epithelial cells donor1 tech_rep1 : CNhs1258": "epithelial cells (nasal)",
    "nasal epithelial cells donor2 : CNhs1257": "epithelial cells (nasal)",
    "Natural Killer Cells donor1 : CNhs1085": "natural killer cells",
    "Natural Killer Cells donor2 : CNhs1195": "natural killer cells",
    "Natural Killer Cells donor3 : CNhs1200": "natural killer cells",
    "Neural stem cells donor1 : CNhs1106": "neural stem cells",
    "Neural stem cells donor2 : CNhs1138": "neural stem cells",
    "neuroblastoma cell line:CHP-134 : CNhs1127": "CHP-134",
    "neuroblastoma cell line:NB-1 : CNhs1128": "NB-1",
    "neuroblastoma cell line:NBsusSR : CNhs1181": "NBsusSR",
    "neuroblastoma cell line:NH-12 : CNhs1181": "NH-12",
    "neuroectodermal tumor cell line:FU-RPNT-1 : CNhs1174": "FU-RPNT-1",
    "neuroectodermal tumor cell line:FU-RPNT-2 : CNhs1175": "FU-RPNT-2",
    "neuroectodermal tumor  cell line:TASK1 : CNhs1186": "TASK1",
    "neuroepithelioma cell line:SK-N-MC : CNhs1185": "SK-N-MC",
    "neurofibroma cell line:Hs 53.T : CNhs1185": "Hs 53.T",
    "Neurons donor1 : CNhs1233": "neurons",
    "Neurons donor2 : CNhs1272": "neurons",
    "Neurons donor3 : CNhs1381": "neurons",
    "Neutrophils donor1 : CNhs1086": "neutrophils",
    "Neutrophils donor2 : CNhs1195": "neutrophils",
    "Neutrophils donor3 : CNhs1190": "neutrophils",
    "NK T cell leukemia cell line:KHYG-1 : CNhs1186": "KHYG-1",
    "non T non B acute lymphoblastic leukemia (ALL) cell line:P30/OHK : CNhs1074": "P30/OHK",
    "normal embryonic palatal mesenchymal cell line:HEPM : CNhs1189": "HEPM",
    "normal intestinal epithelial cell line:FHs 74 Int : CNhs1195": "FHs 74 Int",
    "nucleus accumbens adult pool1 : CNhs1064": "nucleus accumbens",
    "Nucleus Pulposus Cell donor1 : CNhs1088": "nucleus pulposus cells",
    "Nucleus Pulposus Cell donor2 : CNhs1201": "nucleus pulposus cells",
    "Nucleus Pulposus Cell donor3 : CNhs1206": "nucleus pulposus cells",
    "occipital cortex - adult donor10196 : CNhs1379": "occipital cortex",
    "occipital cortex adult donor10252 : CNhs1232": "occipital cortex",
    "occipital lobe adult donor1 : CNhs1178": "occipital lobe",
    "occipital lobe fetal donor1 : CNhs1178": "occipital lobe (fetal)",
    "occipital pole adult pool1 : CNhs1064": "occipital pole",
    "Olfactory epithelial cells donor1 : CNhs1381": "olfactory epithelium",
    "Olfactory epithelial cells donor2 : CNhs1381": "olfactory epithelium",
    "Olfactory epithelial cells donor3 : CNhs1381": "olfactory epithelium",
    "Olfactory epithelial cells donor4 : CNhs1381": "olfactory epithelium",
    "olfactory region adult : CNhs1261": "olfactory region",
    "optic nerve donor1 : CNhs1344": "optic nerve",
    "oral squamous cell carcinoma cell line:Ca9-22 : CNhs1075": "Ca9-22",
    "oral squamous cell carcinoma cell line:HO-1-u-1 : CNhs1128": "HO-1-u-1",
    "oral squamous cell carcinoma cell line:HSC-3 : CNhs1171": "HSC-3",
    "oral squamous cell carcinoma cell line:SAS : CNhs1181": "SAS",
    "Osteoblast - differentiated donor1 : CNhs1131": "osteoblasts (differentiated)",
    "Osteoblast - differentiated donor2 : CNhs1198": "osteoblasts (differentiated)",
    "Osteoblast - differentiated donor3 : CNhs1203": "osteoblasts (differentiated)",
    "Osteoblast donor2 : CNhs1138": "osteoblasts",
    "Osteoblast donor3 : CNhs1203": "osteoblasts",
    "osteoclastoma cell line:Hs 706.T : CNhs1183": "Hs 706.T",
    "osteosarcoma cell line:143B/TK^(-)neo^(R) : CNhs1127": "143B/TK^(-)neo^(R)",
    "osteosarcoma cell line:HS-Os-1 : CNhs1129": "HS-Os-1",
    "ovary adult pool1 : CNhs1062": "ovary",
    "pagetoid sarcoma cell line:Hs 925.T : CNhs1185": "Hs 925.T",
    "pancreas adult donor1 : CNhs1175": "pancreas",
    "pancreatic carcinoma cell line:NOR-P1 : CNhs1183": "NOR-P1",
    "Pancreatic stromal cells donor1 : CNhs1087": "stromal cells (pancreas)",
    "papillary adenocarcinoma cell line:8505C : CNhs1171": "8505C",
    "papillotubular adenocarcinoma cell line:TGBC18TKB : CNhs1073": "TGBC18TKB",
    "paracentral gyrus adult pool1 : CNhs1064": "paracentral gyrus",
    "parietal lobe - adult donor10196 : CNhs1379": "parietal lobe",
    "parietal lobe adult donor10252 : CNhs1231": "parietal lobe",
    "parietal lobe adult pool1 : CNhs1064": "parietal lobe",
    "parietal lobe fetal donor1 : CNhs1178": "parietal lobe (fetal)",
    "parotid gland adult : CNhs1284": "parotid gland",
    "penis adult : CNhs1285": "penis",
    "Pericytes donor1 : CNhs1131": "pericytes",
    "Pericytes donor2 : CNhs1207": "pericytes",
    "Peripheral Blood Mononuclear Cells donor1 : CNhs1086": "peripheral blood mononuclear cells",
    "Peripheral Blood Mononuclear Cells donor2 : CNhs1195": "peripheral blood mononuclear cells",
    "Peripheral Blood Mononuclear Cells donor3 : CNhs1200": "peripheral blood mononuclear cells",
    "peripheral neuroectodermal tumor cell line:KU-SN : CNhs1183": "KU-SN",
    "pharyngeal carcinoma cell line:Detroit 562 : CNhs1184": "Detroit 562",
    "pineal gland - adult donor10196 : CNhs1380": "pineal gland",
    "pineal gland adult donor10252 : CNhs1222": "pineal gland",
    "pituitary gland - adult donor10196 : CNhs1380": "pituitary gland",
    "pituitary gland adult donor10252 : CNhs1222": "pituitary gland",
    "placenta adult pool1 : CNhs1062": "placenta",
    "Placental Epithelial Cells donor1 : CNhs1107": "epithelial cells (placenta)",
    "Placental Epithelial Cells donor2 : CNhs1138": "epithelial cells (placenta)",
    "Placental Epithelial Cells donor3 : CNhs1203": "epithelial cells (placenta)",
    "plasma cell leukemia cell line:ARH-77 : CNhs1280": "ARH-77",
    "pleomorphic hepatocellular carcinoma cell line:SNU-387 : CNhs1193": "SNU-387",
    "pons adult pool1 : CNhs1064": "pons",
    "postcentral gyrus adult pool1 : CNhs1063": "postcentral gyrus",
    "Preadipocyte - breast donor1 : CNhs1105": "preadipocytes (breast)",
    "Preadipocyte - breast donor2 : CNhs1197": "preadipocytes (breast)",
    "Preadipocyte - omental donor1 : CNhs1106": "preadipocytes (omentum)",
    "Preadipocyte - omental donor2 : CNhs1190": "preadipocytes (omentum)",
    "Preadipocyte - omental donor3 : CNhs1201": "preadipocytes (omentum)",
    "Preadipocyte - perirenal donor1 : CNhs1206": "preadipocytes (perirenal)",
    "Preadipocyte - subcutaneous donor2 : CNhs1198": "preadipocytes (subcutaneous)",
    "Preadipocyte - subcutaneous donor3 : CNhs1203": "preadipocytes (subcutaneous)",
    "Preadipocyte - visceral donor1 : CNhs1108": "preadipocytes (viscera)",
    "Preadipocyte - visceral donor2 : CNhs1198": "preadipocytes (viscera)",
    "Preadipocyte - visceral donor3 : CNhs1203": "preadipocytes (viscera)",
    "prostate adult pool1 : CNhs1062": "prostate",
    "prostate cancer cell line:DU145 : CNhs1126": "DU145",
    "prostate cancer cell line:PC-3 : CNhs1124": "PC-3",
    "Prostate Epithelial Cells donor2 : CNhs1197": "epithelial cells (prostate)",
    "Prostate Epithelial Cells donor3 : CNhs1201": "epithelial cells (prostate)",
#    "Prostate Epithelial Cells (polarized) donor1 : CNhs1088": "epithelial cells (prostate)",
    "Prostate Stromal Cells donor1 : CNhs1088": "stromal cells (prostate)",
    "Prostate Stromal Cells donor2 : CNhs1197": "stromal cells (prostate)",
    "Prostate Stromal Cells donor3 : CNhs1201": "stromal cells (prostate)",
    "putamen adult donor10196 : CNhs1232": "putamen",
    "rectal cancer cell line:TT1TKB : CNhs1125": "TT1TKB",
    "rectum fetal donor1 : CNhs1177": "rectum (fetal)",
    "renal cell carcinoma cell line:OS-RC-2 : CNhs1072": "OS-RC-2",
    "renal cell carcinoma cell line:TUHR10TKB : CNhs1125": "TUHR10TKB",
    "Renal Cortical Epithelial Cells donor1 : CNhs1133": "epithelial cells (renal cortex)",
    "Renal Cortical Epithelial Cells donor2 : CNhs1272": "epithelial cells (renal cortex)",
    "Renal Epithelial Cells donor1 : CNhs1133": "epithelial cells (renal)",
    "Renal Epithelial Cells donor2 : CNhs1208": "epithelial cells (renal)",
    "Renal Epithelial Cells donor3 : CNhs1273": "epithelial cells (renal)",
    "Renal Glomerular Endothelial Cells donor1 : CNhs1207": "endothelial cells (renal glomerulus)",
    "Renal Glomerular Endothelial Cells donor2 : CNhs1208": "endothelial cells (renal glomerulus)",
    "Renal Glomerular Endothelial Cells donor3 : CNhs1262": "endothelial cells (renal glomerulus)",
    "Renal Glomerular Endothelial Cells donor4 : CNhs1308": "endothelial cells (renal glomerulus)",
    "Renal Mesangial Cells donor1 : CNhs1133": "mesangial cells",
    "Renal Mesangial Cells donor3 : CNhs1212": "mesangial cells",
    "Renal Proximal Tubular Epithelial Cell donor1 : CNhs1133": "epithelial cell (renal proximal tubule)",
    "Renal Proximal Tubular Epithelial Cell donor2 : CNhs1208": "epithelial cell (renal proximal tubule)",
    "Renal Proximal Tubular Epithelial Cell donor3 : CNhs1212": "epithelial cell (renal proximal tubule)",
    "Reticulocytes biol_ rep1 : CNhs1355": "reticulocytes",
    "Reticulocytes biol_ rep2 : CNhs1355": "reticulocytes",
    "retina adult pool1 : CNhs1063": "retina",
    "Retinal Pigment Epithelial Cells donor0 : CNhs1084": "retinal pigment epithelium",
    "Retinal Pigment Epithelial Cells donor1 : CNhs1133": "retinal pigment epithelium",
    "Retinal Pigment Epithelial Cells donor3 : CNhs1273": "retinal pigment epithelium",
    "retinoblastoma cell line:Y79 : CNhs1126": "Y79",
    "rhabdomyosarcoma cell line:KYM-1 : CNhs1187": "KYM-1",
    "rhabdomyosarcoma cell line:RMS-YM : CNhs1126": "RMS-YM",
#    "SABiosciences XpressRef Human Universal Total RNA pool1 : CNhs1061": ,
    "sacrococcigeal teratoma cell line:HTST : CNhs1182": "HTST",
    "salivary acinar cells donor1 : CNhs1281": "acinar cells (salivary)",
    "salivary acinar cells donor2 : CNhs1281": "acinar cells (salivary)",
    "salivary acinar cells donor3 : CNhs1281": "acinar cells (salivary)",
    "salivary gland adult pool1 : CNhs1167": "salivary gland",
    "schwannoma cell line:HS-PSS : CNhs1118": "HS-PSS",
    "schwannoma cell line:HS-PSS tech_rep2 : CNhs1124": "HS-PSS",
    "Sebocyte donor1 : CNhs1084": "sebocytes",
    "Sebocyte donor2 : CNhs1195": "sebocytes",
    "seminal vesicle adult : CNhs1285": "seminal vesicle",
    "serous adenocarcinoma cell line:JHOS-2 : CNhs1174": "JHOS-2",
#    "serous adenocarcinoma cell line:SK-OV-3-R after co-culture with SOC-57-02-G biol_rep1 : CNhs1350": "SK-OV-3-R",
    "serous adenocarcinoma cell line:SK-OV-3-R biol_rep1 : CNhs1309": "SK-OV-3-R",
    "serous cystadenocarcinoma cell line:HTOA : CNhs1182": "HTOA",
    "Sertoli Cells donor1 : CNhs1085": "Sertoli cells",
    "signet ring carcinoma cell line:Kato III : CNhs1075": "Kato III",
    "signet ring carcinoma cell line:NUGC-4 : CNhs1127": "NUGC-4",
    "skeletal muscle adult pool1 : CNhs1062": "skeletal muscle",
    "Skeletal muscle cells differentiated into Myotubes - multinucleated donor1 : CNhs1108": "myotubes (differentiated from skeletal muscle cells)",
    "Skeletal Muscle Cells donor1 : CNhs1108": "skeletal muscle cells",
    "Skeletal Muscle Cells donor4 : CNhs1205": "skeletal muscle cells",
    "Skeletal Muscle Cells donor5 : CNhs1205": "skeletal muscle cells",
    "Skeletal Muscle Cells donor6 : CNhs1206": "skeletal muscle cells",
    "skeletal muscle fetal donor1 : CNhs1177": "skeletal muscle (fetal)",
    "Skeletal Muscle Satellite Cells donor1 : CNhs1086": "skeletal muscle satellite cells",
    "Skeletal Muscle Satellite Cells donor2 : CNhs1196": "skeletal muscle satellite cells",
    "Skeletal Muscle Satellite Cells donor3 : CNhs1200": "skeletal muscle satellite cells",
    "skeletal muscle - soleus muscle donor1 : CNhs1345": "skeletal muscle (soleus)",
    "skin fetal donor1 : CNhs1177": "skin (fetal)",
    "Small Airway Epithelial Cells donor1 : CNhs1088": "epithelial cells (small airways)",
    "Small Airway Epithelial Cells donor2 : CNhs1197": "epithelial cells (small airways)",
    "Small Airway Epithelial Cells donor3 : CNhs1201": "epithelial cells (small airways)",
    "small cell cervical cancer cell line:HCSC-1 : CNhs1188": "HCSC-1",
    "small cell gastrointestinal carcinoma cell line:ECC10 : CNhs1173": "ECC10",
    "small-cell gastrointestinal carcinoma cell line:ECC4 : CNhs1173": "ECC4",
    "small cell lung carcinoma cell line:DMS 144 : CNhs1280": "DMS 144",
    "small cell lung carcinoma cell line:LK-2 : CNhs1128": "LK-2",
    "small cell lung carcinoma cell line:NCI-H82 : CNhs1280": "NCI-H82",
    "small cell lung carcinoma cell line:WA-hT : CNhs1181": "WA-hT",
    "small intestine adult pool1 : CNhs1063": "small intestine",
    "small intestine fetal donor1 : CNhs1177": "small intestine (fetal)",
    "smooth muscle adult pool1 : CNhs1175": "smooth muscle",
    "Smooth Muscle Cells - Aortic donor0 : CNhs1083": "smooth muscle cells (aorta)",
    "Smooth Muscle Cells - Aortic donor1 : CNhs1108": "smooth muscle cells (aorta)",
    "Smooth Muscle Cells - Aortic donor2 : CNhs1130": "smooth muscle cells (aorta)",
    "Smooth Muscle Cells - Aortic donor3 : CNhs1130": "smooth muscle cells (aorta)",
    "Smooth Muscle Cells - Brachiocephalic donor1 : CNhs1108": "smooth muscle cells (brachiocephalic artery)",
    "Smooth Muscle Cells - Brachiocephalic donor3 : CNhs1204": "smooth muscle cells (brachiocephalic artery)",
    "Smooth Muscle Cells - Brain Vascular donor1 : CNhs1086": "smooth muscle cells (brain vasculature)",
    "Smooth Muscle Cells - Brain Vascular donor2 : CNhs1190": "smooth muscle cells (brain vasculature)",
    "Smooth Muscle Cells - Brain Vascular donor3 : CNhs1200": "smooth muscle cells (brain vasculature)",
    "Smooth Muscle Cells - Bronchial donor1 : CNhs1132": "smooth muscle cells (bronchi)",
    "Smooth Muscle Cells - Bronchial donor2 : CNhs1234": "smooth muscle cells (bronchi)",
    "Smooth Muscle Cells - Carotid donor1 : CNhs1108": "smooth muscle cells (carotid)",
    "Smooth Muscle Cells - Carotid donor3 : CNhs1204": "smooth muscle cells (carotid)",
    "Smooth Muscle Cells - Colonic donor1 : CNhs1086": "smooth muscle cells (colon)",
    "Smooth Muscle Cells - Colonic donor2 : CNhs1196": "smooth muscle cells (colon)",
    "Smooth Muscle Cells - Colonic donor3 : CNhs1200": "smooth muscle cells (colon)",
    "Smooth Muscle Cells - Coronary Artery donor1 : CNhs1108": "smooth muscle cells (coronary artery)",
    "Smooth Muscle Cells - Coronary Artery donor2 : CNhs1198": "smooth muscle cells (coronary artery)",
    "Smooth Muscle Cells - Coronary Artery donor3 : CNhs1204": "smooth muscle cells (coronary artery)",
    "Smooth Muscle Cells - Esophageal donor1 : CNhs1132": "smooth muscle cells (esophagus)",
    "Smooth Muscle Cells - Internal Thoracic Artery donor2 : CNhs1198": "smooth muscle cells (internal thoracic artery)",
    "Smooth Muscle Cells - Internal Thoracic Artery donor3 : CNhs1204": "smooth muscle cells (internal thoracic artery)",
    "Smooth Muscle Cells - Prostate donor1 : CNhs1192": "smooth muscle cells (prostate)",
    "Smooth Muscle Cells - Prostate donor2 : CNhs1197": "smooth muscle cells (prostate)",
    "Smooth Muscle Cells - Pulmonary Artery donor2 : CNhs1198": "smooth muscle cells (pulmonary artery)",
    "Smooth Muscle Cells - Subclavian Artery donor1 : CNhs1109": "smooth muscle cells (subclavian artery)",
    "Smooth Muscle Cells - Subclavian Artery donor2 : CNhs1199": "smooth muscle cells (subclavian artery)",
    "Smooth Muscle Cells - Subclavian Artery donor3 : CNhs1204": "smooth muscle cells (subclavian artery)",
    "Smooth Muscle Cells - Tracheal donor1 : CNhs1132": "smooth muscle cells (trachea)",
    "Smooth Muscle Cells - Tracheal donor3 : CNhs1289": "smooth muscle cells (trachea)",
    "Smooth Muscle Cells - Umbilical artery donor0 : CNhs1083": "smooth muscle cells (umbilical artery)",
    "Smooth Muscle Cells - Umbilical Artery donor1 : CNhs1109": "smooth muscle cells (umbilical artery)",
    "Smooth Muscle Cells - Umbilical Artery donor2 : CNhs1199": "smooth muscle cells (umbilical artery)",
    "Smooth Muscle Cells - Umbilical Artery donor3 : CNhs1204": "smooth muscle cells (umbilical artery)",
    "Smooth Muscle Cells - Umbilical Vein donor1 : CNhs1259": "smooth muscle cells (umbilical vein)",
    "Smooth Muscle Cells - Umbilical Vein donor2 : CNhs1256": "smooth muscle cells (umbilical vein)",
    "Smooth Muscle Cells - Uterine donor3 : CNhs1192": "smooth muscle cells (uterus)",
    "somatostatinoma cell line:QGP-1 : CNhs1186": "QGP-1",
    "spinal cord - adult donor10196 : CNhs1380": "spinal cord",
    "spinal cord adult donor10252 : CNhs1222": "spinal cord",
    "spinal cord fetal donor1 : CNhs1176": "spinal cord (fetal)",
    "spindle cell sarcoma cell line:Hs 132.T : CNhs1185": "Hs 132.T",
    "spleen adult pool1 : CNhs1063": "spleen",
    "spleen fetal pool1 : CNhs1065": "spleen (fetal)",
    "splenic lymphoma with villous lymphocytes cell line:SLVL : CNhs1074": "SLVL",
    "squamous cell carcinoma cell line:EC-GI-10 : CNhs1125": "EC-GI-10",
    "squamous cell carcinoma cell line:T3M-5 : CNhs1173": "T3M-5",
    "squamous cell lung carcinoma cell line:EBC-1 : CNhs1127": "EBC-1",
    "stomach fetal donor1 : CNhs1177": "stomach (fetal)",
    "submaxillary gland adult : CNhs1285": "submaxillary gland",
    "substantia nigra adult donor10252 : CNhs1231": "substantia nigra",
    "synovial sarcoma cell line:HS-SY-II : CNhs1124": "HS-SY-II",
    "Synoviocyte donor2 : CNhs1199": "synoviocytes",
    "Synoviocyte donor3 : CNhs1205": "synoviocytes",
    "temporal lobe adult pool1 : CNhs1063": "temporal lobe",
    "temporal lobe fetal donor1 tech_rep1 : CNhs1177": "temporal lobe (fetal)",
    "temporal lobe fetal donor1 tech_rep2 : CNhs1299": "temporal lobe (fetal)",
    "tenocyte donor1 : CNhs1263": "tenocytes",
    "tenocyte donor2 : CNhs1264": "tenocytes",
    "tenocyte donor3 : CNhs1264": "tenocytes",
    "teratocarcinoma cell line:NCC-IT-A3 : CNhs1187": "NCC-IT-A3",
    "teratocarcinoma cell line:NCR-G1 : CNhs1188": "NCR-G1",
    "teratocarcinoma  cell line:PA-1 : CNhs1189": "PA-1",
    "testicular germ cell embryonal carcinoma cell line:ITO-II : CNhs1187": "ITO-II",
    "testicular germ cell embryonal carcinoma cell line:NEC14 : CNhs1235": "NEC14",
    "testicular germ cell embryonal carcinoma cell line:NEC15 : CNhs1236": "NEC15",
    "testicular germ cell embryonal carcinoma cell line:NEC8 : CNhs1172": "NEC8",
    "testis adult pool1 : CNhs1063": "testis",
    "testis adult pool2 : CNhs1299": "testis",
    "thalamus - adult donor10196 : CNhs1379": "thalamus",
    "thalamus adult donor10252 : CNhs1231": "thalamus",
    "throat adult : CNhs1285": "throat",
    "throat fetal donor1 : CNhs1177": "throat (fetal)",
    "thymus adult pool1 : CNhs1063": "thymus",
    "thymus fetal pool1 : CNhs1065": "thymus (fetal)",
    "thyroid adult pool1 : CNhs1063": "thyroid",
    "thyroid carcinoma cell line:TCO-1 : CNhs1187": "TCO-1",
    "thyroid fetal donor1 : CNhs1176": "thyroid (fetal)",
    "tongue adult : CNhs1285": "tongue",
    "tongue fetal donor1 : CNhs1176": "tongue (fetal)",
    "tonsil adult pool1 : CNhs1065": "tonsils",
    "Trabecular Meshwork Cells donor1 : CNhs1134": "trabecular meshwork cells",
    "Trabecular Meshwork Cells donor3 : CNhs1212": "trabecular meshwork cells",
    "trachea adult pool1 : CNhs1063": "trachea",
    "trachea fetal donor1 : CNhs1176": "trachea (fetal)",
    "Tracheal Epithelial Cells donor1 : CNhs1109": "epithelial cells (trachea)",
    "Tracheal Epithelial Cells donor2 : CNhs1199": "epithelial cells (trachea)",
    "Tracheal Epithelial Cells donor3 : CNhs1205": "epithelial cells (trachea)",
    "transitional-cell carcinoma cell line:5637 : CNhs1073": "5637",
    "transitional-cell carcinoma cell line:JMSU1 : CNhs1126": "JMSU1",
    "tridermal teratoma cell line:HGRT : CNhs1182": "HGRT",
    "tubular adenocarcinoma cell line:SUIT-2 : CNhs1188": "SUIT-2",
    "umbilical cord fetal donor1 : CNhs1176": "umbilical cord (fetal)",
#    "Universal RNA - Human Normal Tissues Biochain pool1 : CNhs1061": ,
    "Urothelial cells donor0 : CNhs1084": "urothelial cells",
    "Urothelial Cells donor1 : CNhs1133": "urothelial cells",
    "Urothelial Cells donor2 : CNhs1209": "urothelial cells",
    "Urothelial Cells donor3 : CNhs1212": "urothelial cells",
    "uterus adult pool1 : CNhs1167": "uterus",
    "uterus fetal donor1 : CNhs1176": "uterus (fetal)",
    "vagina adult : CNhs1285": "vagina",
    "vein adult : CNhs1284": "vein",
#    "Whole blood (ribopure) donor090309 donation1 : CNhs1167": "blood (ribopure)",
#    "Whole blood (ribopure) donor090309 donation2 : CNhs1167": "blood (ribopure)",
#    "Whole blood (ribopure) donor090309 donation3 : CNhs1194": "blood (ribopure)",
#    "Whole blood (ribopure) donor090325 donation1 : CNhs1107": "blood (ribopure)",
#    "Whole blood (ribopure) donor090325 donation2 : CNhs1107": "blood (ribopure)",
#    "Whole blood (ribopure) donor090612 donation1 : CNhs1167": "blood (ribopure)",
#    "Whole blood (ribopure) donor090612 donation2 : CNhs1167": "blood (ribopure)",
#    "Whole blood (ribopure) donor090612 donation3 : CNhs1194": "blood (ribopure)",
    "Wilms' tumor cell line:G-401 : CNhs1189": "G-401",
    "Wilms' tumor cell line:HFWT : CNhs1172": "HFWT",
    "xeroderma pigentosum b cell line:XPL 17 : CNhs1181": "XPL 17",
}

grouped_tss_sample_names = {
    
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

    parser = argparse.ArgumentParser(description="this script inserts \"enhancer\" or \"tss\" data from FANATOM into GUD. download \"http://slidebase.binf.ku.dk/human_enhancers/presets/serve/hg19_permissive_enhancers_expression_rle_tpm.csv.gz\" and \"\" for enhancer and tss data, respectively.")

    parser.add_argument("matrix", help="Expression (TPM/RLE normalized) matrix across all FANTOM libraries")

    feats = ["enhancer", "tss"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature", metavar="feature_type")

    # Optional args
    parser.add_argument("-b", "--bed", help="BED file of features on which to focus (e.g. \"robust_enhancers.bed\")")
    parser.add_argument("-g", "--group", action="store_true", help="Group FANTOM samples (default = False)")
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
    feat_type, source_name, bed_file=None, group=False):

    # Initialize
    original_sample_names = []
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
        if not engine.has_table("enhancer"):
            raise ValueError("GUD db does not have \"enhancer\" table!")
        table = Enhancer()
    if feat_type == "tss":
        if not engine.has_table("tss"):
            raise ValueError("GUD db does not have \"tss\" table!")
        table = TSS()
    table.metadata.bind = engine
    table.metadata.create_all(engine)
    mapper(Model, table.__table__)

    # If BED file...
    if bed_file:
        # If BED file exists...
        if os.path.exists(bed_file):
            try:
                # Create BED object
                bed_obj = pybedtools.BedTool(bed_file)
            except:
                warnings.warn("\nCould not read file: \"%s\"\n\tSkipping file...\n" % file_name)

    
    if feat_type == "enhancer":
        # For each line...
        for line in GUDglobals.parse_csv_file(matrix_file, gz):
            # If no samples...
            if len(original_sample_names) == 0:
                for sample in line[1:]:
                    original_sample_names.append(sample[1:-2])
            # ... Else...
            else:
                # Initialize
                samples = {}
                # Get chrom, start, end
                m = re.search("(chr\S+)\:(\d+)\-(\d+)", line.pop(0))
                chrom = m.group(1)
                start = int(m.group(2))
                end = int(m.group(3))
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
                # For each sample...
                for i in range(len(line)):
                    # Initialize
                    cages = float(line[i])
                    sample = original_sample_names[i]
                    # Group samples
                    if group:
                        if sample in grouped_sample_names:
                            samples.setdefault(grouped_sample_names[sample], [])
                            samples[grouped_sample_names[sample]].append(cages)
                    else: samples.setdefault(sample, [cages])
                # For each sample...
                for sample in samples:
                    model.cell_or_tissue = sample
                    # Skip enhancers with 0 cages
                    if sum(samples[sample]) == 0:
                        continue
                    # Upsert model & commit
                    session.merge(model)
                    session.commit()

    if feat_type == "tss":
        # For each line...
        for line in GUDglobals.parse_tsv_file(matrix_file, gz):
            # Skip comments
            if line[0].startswith("#"): continue
            # If no samples...
            if len(original_sample_names) == 0:
                for sample in line[1:]:
                    original_sample_names.append(unquote(sample))
            # ... Else...
            else: pass
            for sample in original_sample_names:
                print(sample)
                m = re.search("(CNhs\d+)", sample)
                for curated_sample in grouped_sample_names:
                    print(curated_sample)
                    n = re.search("(CNhs\d+)", curated_sample)
                    if m.group(1) == n.group(1):
                        print("\"\": \"\"," % (sample, curated_sample))
                        break
            exit(0)
#                # Initialize
#                samples = {}
#                # Get chrom, start, end
#                m = re.search("(chr\S+)\:(\d+)\-(\d+)", line.pop(0))
#                chrom = m.group(1)
#                start = int(m.group(2))
#                end = int(m.group(3))
#                # Ignore non-standard chroms, scaffolds, etc.
#                m = re.search("^chr(\w{1,2})$", chrom)
#                if not m.group(1) in GUDglobals.chroms: continue
#                # Create model
#                model = Model()
#                model.bin = assign_bin(int(start), int(end))
#                model.chrom = chrom
#                model.start = start
#                model.end = end
#                model.experiment_type = "CAGE"
#                model.source_name = source_name
#                model.date = today
#                # For each sample...
#                for i in range(len(line)):
#                    # Initialize
#                    cages = float(line[i])
#                    sample = original_sample_names[i]
#                    # Group samples
#                    if group:
#                        if sample in grouped_sample_names:
#                            samples.setdefault(grouped_sample_names[sample], [])
#                            samples[grouped_sample_names[sample]].append(cages)
#                    else: samples.setdefault(sample, [cages])
#                # For each sample...
#                for sample in samples:
#                    model.cell_or_tissue = sample
#                    # Skip enhancers with 0 cages
#                    if sum(samples[sample]) == 0:
#                        continue
#                    # Upsert model & commit
#                    session.merge(model)
#                    session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert FANTOM data to GUD database
    insert_fantom_to_gud_db(args.user, args.host, args.port,
        args.db, args.matrix, args.feat_type, args.source,
        args.bed, args.group)