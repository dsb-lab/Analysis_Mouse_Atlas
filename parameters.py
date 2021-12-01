#All the human chosen parameters of the analysis

# QC measures 
#!!! None of this parameters are really chosen in our analysis and are just sloppy bounds of the QC in Pijuan
MIN_COUNTS = 5000
MAX_COUNTS = 10E8
MIN_GENES = 0
MAX_GENES = 10E8
MIN_CELLS = 1
MAX_CELLS = 10E8
MIN_MT_FRACTION = 0
MAX_MT_FRACTION = 1

# QC Doublets and Stripped cells
DOUBLET_LOUVAIN_RESOLUTION = 10  
##Doublets
MAX_N_COUNTS = 39000
DOUBLET_SCORE_MAX_THRESHOLD = 0.15
##Stripped
STRIPPED = 0.0022

#General parameters
## HVGs
HVG_METHOD = "seurat"
N_HVG = None
## PCs
USE_HVGs = True
N_PCS = 50
## Neighbours
METRIC = "correlation"
N_NEIGBOURS = 25
## Louvain
LOUVAIN_INITIALISATION = 2 #Used for initialising the UMAPS

#SAMPLE PARAMETERS
##HVGs
RECOMPUTE_HVGS = False
##PCs
RECOMPUTE_PCS = True
##Cluster threshold
CLUSTER_THRESHOLD = 0.25
##Remove outliers
DISCONNECTION_DISTANCE = 7
MIN_CLUSTER_SIZE = 15

#EXCLUDE FROM HVGs GENES
SKIP_GENES_SEX = ["Xist"]
SKIP_GENES_S_PHASE = \
    ['Mcm5', 'Pcna', 'Tyms', 'Fen1', 'Mcm2', 'Mcm4', 'Rrm1', 'Ung', 'Gins2',
     'Mcm6', 'Cdca7', 'Dtl', 'Prim1', 'Uhrf1', 'Mlf1ip', 'Hells', 'Rfc2',
     'Rpa2', 'Nasp', 'Rad51ap1', 'Gmnn', 'Wdr76', 'Slbp', 'Ccne2', 'Ubr7',
     'Pold3', 'Msh2', 'Atad2', 'Rad51', 'Rrm2', 'Cdc45', 'Cdc6', 'Exo1', 'Tipin',
     'Dscc1', 'Blm', 'Casp8ap2', 'Usp1', 'Clspn', 'Pola1', 'Chaf1b', 'Brip1', 'E2f8']
SKIP_GENES_G2M = \
    ['Hmgb2', 'Cdk1', 'Nusap1', 'Ube2c', 'Birc5', 'Tpx2', 'Top2a', 'Ndc80',
     'Cks2', 'Nuf2', 'Cks1b', 'Mki67', 'Tmpo', 'Cenpf', 'Tacc3', 'Fam64a',
     'Smc4', 'Ccnb2', 'Ckap2l', 'Ckap2', 'Aurkb', 'Bub1', 'Kif11', 'Anp32e',
     'Tubb4b', 'Gtse1', 'Kif20b', 'Hjurp', 'Cdca3', 'Hn1', 'Cdc20', 'Ttk',
     'Cdc25c', 'Kif2c', 'Rangap1', 'Ncapd2', 'Dlgap5', 'Cdca2', 'Cdca8',
     'Ect2', 'Kif23', 'Hmmr', 'Aurka', 'Psrc1', 'Anln', 'Lbr', 'Ckap5',
     'Cenpe', 'Ctcf', 'Nek2', 'G2e3', 'Gas2l3', 'Cbx5', 'Cenpa']
SKIP_GENES = SKIP_GENES_SEX + SKIP_GENES_S_PHASE + SKIP_GENES_G2M

#ANNOTATION
#Louvain resolution
LOUVAIN_RESOLUTION = 1
DPC_CUT_PARAMETERS = [20,0.2]
#Genes used for manual annotation
geneList = [
 'Adh1a2',
 'Bmp4','Bmp8',
 'Cdh1','Cdh2','Cdh11','Cer1',
 'Dlx5','Dmrt2',
 'Ebf1','Egr2','En1','Eomes','Epcam','Eya1','Eya2','Evx1',
 'Fgf5','Foxa2','Foxc1','Foxd1','Foxf1','Foxf2','Foxp2',
 'Hand1','Hand2','Hesx1','Hoxa1','Hoxa2','Hoxa3','Hoxa4','Hoxa10',
 'Irx3','Isl1',
 'Kdr',
 'Lefty1','Lef1','Lfng','Lhx1',
 'Meox1','Mesp1','Mesp2','Msc','Msgn','Myf5',
 'Nanog','Nkx2-5','Noto',
 'Olig2','Osr1','Otx2',
 'Pax1','Pax3','Pax6','Pitx2','Prdm1','Pou3f1','Pou5f1',
 'Ripply2',
 'Shh','Six1','Six3','Snai1','Snai2','Sox1','Sox2','Sox3','Sox9','Sox17',
 'T','Tfap2a','Tfap2c','Tbx1','Tbx4','Tbx5','Tbx6','Tbx18','Tcf15','Tcf21','Twist1','Twist2',
 'Uncx',
 'Wnt1','Wnt3','Wnt3a'
]

#COMPARATIVE GENES
genes_comparative = ["Osr1","Tbx6","Foxf1"]
stages_comparative = {"E7.5":["3","4","7","10","11"],
                      "E8.0":["4","5","6"],
                      "E8.5":["1","2"]}