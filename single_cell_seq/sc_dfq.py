"""
Description: some test.
Date: 2019年12月2日13:21:21
Author: dfq
"""

import os
import scanpy as sc
import numpy as np
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt

# import scipy as sp
from matplotlib import rcParams
from matplotlib import colors
# from gprofiler import gprofiler
#
# import rpy2.rinterface_lib.callbacks
# import logging
# from rpy2.robjects import pandas2ri
# import anndata2ri

# 导入R包
"""
from rpy2.robjects.packages import importr
scran = importr("scran")
RColorBrewer = importr("RColorBrewer")
slingshot = importr("slingshot")
monocle = importr("monocle")
gam = importr("gam")
clusterExperiment = importr("clusterExperiment")
ggplot2 = importr("ggplot2")
plyr = importr("plyr")
MAST = importr("MAST")

# Ignore R warning messages
#Note: this can be commented out to get more verbose R output
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
# %load_ext rpy2.ipython

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
#sc.set_figure_params(dpi=200, dpi_save=300)
# sc.logging.print_versions()
"""


# Set up data loading
def get_raw_data():
    # Data files  设置文件的识别类型
    sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
    sample_id_strings = ['3', '4', '5', '6', '7', '8']
    # file_base = '/home/fdong/data/sc_data/GSM283657'
    file_base = 'E:\sc_data\GSE92332\GSE92332_RAW\GSM283657\GSM283657'
    exp_string = '_Regional_'
    data_file_end = '_matrix.mtx'
    barcode_file_end = '_barcodes.tsv'
    gene_file_end = '_genes.tsv'
    # cc_genes_file = '~/data/Macosko_cell_cycle_genes.txt'
    cc_genes_file = 'E:\\sc_data\\Macosko_cell_cycle_genes.txt'

    # First data set load & annotation
    #Parse Filenames
    sample = sample_strings.pop(0)
    sample_id = sample_id_strings.pop(0)
    data_file = file_base+sample_id+exp_string+sample+data_file_end
    barcode_file = file_base+sample_id+exp_string+sample+barcode_file_end
    gene_file = file_base+sample_id+exp_string+sample+gene_file_end


    #Load data
    adata = sc.read(data_file, cache=True)    # 主要数据存储对象
    adata = adata.transpose()
    adata.X = adata.X.toarray()


    barcodes = pd.read_csv(barcode_file, header=None, sep='\t')
    genes = pd.read_csv(gene_file, header=None, sep='\t')

    #Annotate data
    barcodes.rename(columns={0:'barcode'}, inplace=True)
    barcodes.set_index('barcode', inplace=True)
    adata.obs = barcodes
    adata.obs['sample'] = [sample]*adata.n_obs
    adata.obs['region'] = [sample.split("_")[0]]*adata.n_obs
    adata.obs['donor'] = [sample.split("_")[1]]*adata.n_obs

    genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
    genes.set_index('gene_symbol', inplace=True)
    adata.var = genes

    for i in range(len(sample_strings)):
        # Parse Filenames
        sample = sample_strings[i]
        sample_id = sample_id_strings[i]
        data_file = file_base + sample_id + exp_string + sample + data_file_end
        barcode_file = file_base + sample_id + exp_string + sample + barcode_file_end
        gene_file = file_base + sample_id + exp_string + sample + gene_file_end

        # Load data
        adata_tmp = sc.read(data_file, cache=True)
        adata_tmp = adata_tmp.transpose()
        adata_tmp.X = adata_tmp.X.toarray()

        barcodes_tmp = pd.read_csv(barcode_file, header=None, sep='\t')
        genes_tmp = pd.read_csv(gene_file, header=None, sep='\t')

        # Annotate data
        barcodes_tmp.rename(columns={0: 'barcode'}, inplace=True)
        barcodes_tmp.set_index('barcode', inplace=True)
        adata_tmp.obs = barcodes_tmp
        adata_tmp.obs['sample'] = [sample] * adata_tmp.n_obs
        adata_tmp.obs['region'] = [sample.split("_")[0]] * adata_tmp.n_obs
        adata_tmp.obs['donor'] = [sample.split("_")[1]] * adata_tmp.n_obs

        genes_tmp.rename(columns={0: 'gene_id', 1: 'gene_symbol'}, inplace=True)
        genes_tmp.set_index('gene_symbol', inplace=True)
        adata_tmp.var = genes_tmp
        adata_tmp.var_names_make_unique()

        # Concatenate to main adata object
        adata = adata.concatenate(adata_tmp, batch_key='sample_id')
        # adata.var['gene_id'] = adata.var['gene_id-1']
        # adata.var.drop(columns=['gene_id-1', 'gene_id-0'], inplace=True)
        adata.obs.drop(columns=['sample_id'], inplace=True)
        adata.obs_names = [c.split("-")[0] for c in adata.obs_names]
        adata.obs_names_make_unique(join='_')


    #Assign variable names and gene id columns
    adata.var_names = [g.split("_")[1] for g in adata.var_names]
    adata.var['gene_id'] = [g.split("_")[1] for g in adata.var['gene_id']]

    return adata

def pre_qc(adata):
    # Quality control - calculate QC covariates
    adata.obs['n_counts'] = adata.X.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes'] = (adata.X > 0).sum(1)

    mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1) / adata.obs['n_counts']  # 计算每个细胞中线粒体基因比例
    return adata

def sc_qc_violin_plot(adata, keys, group='sample', **kwds):
    if keys in adata.obs.columns:
        sc.pl.violin(adata, keys, groupby=group, cut=0, save="_{0}.png".format(keys), show=False, **kwds)
        sc.pl.violin(adata, keys, groupby=group, cut=0, save="_{0}.pdf".format(keys), show=False, **kwds)
    else:
        print("adata.obs.columns没有该{0}数据列".format(keys))

def sc_qc_scatter_plot(adata,x,y,color=None,save="00",**kwds):
    if (x in adata.obs.columns) and (y in adata.obs.columns):
        sc.pl.scatter(adata, x, y, color=color, save="_{0}.png".format(save), show=False, **kwds)
        sc.pl.scatter(adata, x, y, color=color, save="_{0}.pdf".format(save), show=False, **kwds)
    else:
        print("adata.obs.columns没有该{0}or{1}数据列".format(x,y))

def sc_qc_dist_plot(data,keys,*arg,**kwargs):
    fig,ax = plt.subplots()
    sb.distplot(data,*arg,**kwargs)
    fig.savefig("figures/{}.png".format(keys))
    fig.savefig("figures/{}.pdf".format(keys))
    print("figures/{} has been ploted!".format(keys))
    plt.close()

def sc_top_gene_plot(adata,sample_strings):

    def _top_plot(data,name):
        """
        ploting.
        :param data: pd.DataFrame
        :param name: save filename
        :return:
        """
        fig, ax = plt.subplots()
        ax.boxplot(data, vert=False, showfliers=True, labels=hh.index,
                   flierprops={"markersize": 2, "markerfacecolor": "red", "markeredgecolor": "red"})
        ax.set_title(name)
        ax.set_xlabel("percent(%) of total counts")
        fig.savefig("figures/{0}_top20.png".format(name))
        fig.savefig("figures/{0}_top20.pdf".format(name))
        print("{0}_top20  has been saved!".format(name))
        plt.close()

    # sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
    # 生成所有数据集
    aa = pd.DataFrame(adata.X)
    aa.index = adata.obs.index
    aa.columns = adata.var_names

    # select top20 genes for each sample
    for sample in sample_strings:
        # 生成样本的数据集 如：Duo_M2
        dd = aa[adata.obs["sample"] == sample].T
        dd = dd/dd.sum(0)*100
        hh = dd.loc[dd.T.quantile([0.5]).T.sort_values(by=0.5)[-20:].index, :]  # 基因counts比例中位数最高
        _top_plot(hh, sample)

    # select top20 genes for all sample
    hh = aa.T.loc[aa.quantile([0.5]).T.sort_values(by=0.5)[-20:].index, :]
    _top_plot(hh, "all")

def sc_qc_filter(adata, **kwargs):
    before_cells = adata.shape[0]
    print("Total number of cells: [{0}]".format(before_cells))
    # Filter cells according to identified QC thresholds:
    for key,value in kwargs.items():
        if key == "min_counts":
            adata = adata[adata.obs['n_counts'] > value]
        elif key == "min_genes":
            adata = adata[adata.obs['n_genes'] > value]
        elif key == "max_counts":
            adata = adata[adata.obs['n_counts'] < value]
        elif key == "max_genes":
            adata = adata[adata.obs['n_genes'] < value]
        elif key =="mt_frac":
            adata = adata[adata.obs['mt_frac'] < value]
        else:
            print("Have nothing to do with {}={}.".format(key,value))
        after_cells = adata.shape[0]
        gap = before_cells - after_cells
        before_cells = after_cells

        print("Number of cells after {0} filter:[{1}]\tand [{2}] cells filtered out"
              .format(key,after_cells,gap))
    # sc.pp.filter_cells(adata, max_counts = 40000)
    # adata = adata[adata.obs['n_counts'] < 40000]

    sc.pp.filter_genes(adata, min_cells=20)
    return adata

def get_size_factors(adata_pp):

    # 保存数据 供R分析
    adata_pp.obs.to_csv("obs.txt",sep="\t",index_label="cell")
    aa = pd.DataFrame(adata_pp.X)
    # aa.index = adata_pp.obs.index
    aa.columns = adata_pp.var_names
    aa.to_csv("X.txt",sep="\t")
    os.system("Rscript Normalization.R")

    size_factors= pd.read_csv("size_factors.txt",header=0)
    return size_factors

def to_visualization(adata):
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata)

    sc.tl.tsne(adata, n_jobs=12)  # Note n_jobs works for MulticoreTSNE, but not regular implementation)
    sc.tl.umap(adata)
    sc.tl.diffmap(adata)
    sc.tl.draw_graph(adata)

    sc.pl.pca_scatter(adata, color='n_counts',show=False, save="pca.png")
    sc.pl.tsne(adata, color='n_counts',show=False, save="tsne.png")
    sc.pl.umap(adata, color='n_counts',show=False, save="umap.png")
    sc.pl.diffmap(adata, color='n_counts', components=['1,2', '1,3'],show=False, save="diffmap.png")
    sc.pl.draw_graph(adata, color='n_counts',show=False, save="draw_graph.png")

def cell_cycle_scoring(adata,cc_genes_file):
    # Score cell cycle and visualize the effect:
    cc_genes = pd.read_table(cc_genes_file, delimiter='\t')
    s_genes = cc_genes['S'].dropna()
    g2m_genes = cc_genes['G2.M'].dropna()

    s_genes_mm = [gene.lower().capitalize() for gene in s_genes]
    g2m_genes_mm = [gene.lower().capitalize() for gene in g2m_genes]

    s_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
    g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)

    sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False,show=False,save="S_score_G2M_score.png")
    sc.pl.umap(adata, color='phase', use_raw=False,show=False,save="cell_cycle_phase.png")

def to_clustering(adata):
    # Perform clustering - using highly variable genes
    sc.tl.louvain(adata, key_added='louvain_r1')
    sc.tl.louvain(adata, resolution=0.5, key_added='louvain_r0.5', random_state=10)
    adata.obs['louvain_r0.5'].value_counts()

    sc.pl.umap(adata, color=['louvain_r1', 'louvain_r0.5'], palette=sc.pl.palettes.default_64,show=False,save="louvain_r1_r0.5.png")
    sc.pl.umap(adata, color=['region', 'n_counts'],show=False,save="region_n_counts.png")
    sc.pl.umap(adata, color=['log_counts', 'mt_frac'],show=False,save="log_counts_mt_frac.png")

def mark_geenes(adata):
    # Calculate marker genes
    sc.tl.rank_genes_groups(adata, groupby='louvain_r0.5', key_added='rank_genes_r0.5')

    # Plot marker genes
    sc.pl.rank_genes_groups(adata, key='rank_genes_r0.5', groups=['0', '1', '2'], fontsize=12)
    sc.pl.rank_genes_groups(adata, key='rank_genes_r0.5', groups=['3', '4', '5'], fontsize=12)
    sc.pl.rank_genes_groups(adata, key='rank_genes_r0.5', groups=['6', '7', '8'], fontsize=12)


    #Known marker genes:
    marker_genes = dict()
    marker_genes['Stem'] = ['Lgr5', 'Ascl2', 'Slc12a2', 'Axin2', 'Olfm4', 'Gkn3']
    marker_genes['Enterocyte (Proximal)'] = ['Gsta1','Rbp2','Adh6a','Apoa4','Reg3a','Creb3l3','Cyp3a13','Cyp2d26','Ms4a10','Ace','Aldh1a1','Rdh7','H2-Q2', 'Hsd17b6','Gstm3','Gda','Apoc3','Gpd1','Fabp1','Slc5a1','Mme','Cox7a1','Gsta4','Lct','Khk','Mttp','Xdh','Sult1b1', 'Treh','Lpgat1','Dhrs1','Cyp2c66','Ephx2','Cyp2c65','Cyp3a25','Slc2a2','Ugdh','Gstm6','Retsat','Ppap2a','Acsl5', 'Cyb5r3','Cyb5b','Ckmt1','Aldob','Ckb','Scp2','Prap1']
    marker_genes['Enterocyte (Distal)'] = ['Tmigd1','Fabp6','Slc51b','Slc51a','Mep1a','Fam151a','Naaladl1','Slc34a2','Plb1','Nudt4','Dpep1','Pmp22','Xpnpep2','Muc3','Neu1','Clec2h','Phgr1','2200002D01Rik','Prss30','Cubn','Plec','Fgf15','Crip1','Krt20','Dhcr24','Myo15b','Amn','Enpep','Anpep','Slc7a9','Ocm','Anxa2','Aoc1','Ceacam20','Arf6','Abcb1a','Xpnpep1','Vnn1','Cndp2','Nostrin','Slc13a1','Aspa','Maf','Myh14']
    marker_genes['Goblet'] = ['Agr2', 'Fcgbp', 'Tff3', 'Clca1', 'Zg16', 'Tpsg1', 'Muc2', 'Galnt12', 'Atoh1', 'Rep15', 'S100a6', 'Pdia5', 'Klk1', 'Pla2g10', 'Spdef', 'Lrrc26', 'Ccl9', 'Bace2', 'Bcas1', 'Slc12a8', 'Smim14', 'Tspan13', 'Txndc5', 'Creb3l4', 'C1galt1c1', 'Creb3l1', 'Qsox1', 'Guca2a', 'Scin', 'Ern2', 'AW112010', 'Fkbp11', 'Capn9', 'Stard3nl', 'Slc50a1', 'Sdf2l1', 'Hgfa', 'Galnt7', 'Hpd', 'Ttc39a', 'Tmed3', 'Pdia6', 'Uap1', 'Gcnt3', 'Tnfaip8', 'Dnajc10', 'Ergic1', 'Tsta3', 'Kdelr3', 'Foxa3', 'Tpd52', 'Tmed9', 'Spink4', 'Nans', 'Cmtm7', 'Creld2', 'Tm9sf3', 'Wars', 'Smim6', 'Manf', 'Oit1', 'Tram1', 'Kdelr2', 'Xbp1', 'Serp1', 'Vimp', 'Guk1', 'Sh3bgrl3', 'Cmpk1', 'Tmsb10', 'Dap', 'Ostc', 'Ssr4', 'Sec61b', 'Pdia3', 'Gale', 'Klf4', 'Krtcap2', 'Arf4', 'Sep15', 'Ssr2', 'Ramp1', 'Calr', 'Ddost']
    marker_genes['Paneth'] = ['Gm15284', 'AY761184', 'Defa17', 'Gm14851', 'Defa22', 'Defa-rs1', 'Defa3', 'Defa24', 'Defa26', 'Defa21', 'Lyz1', 'Gm15292', 'Mptx2', 'Ang4']
    marker_genes['Enteroendocrine'] = ['Chgb', 'Gfra3', 'Cck', 'Vwa5b2', 'Neurod1', 'Fev', 'Aplp1', 'Scgn', 'Neurog3', 'Resp18', 'Trp53i11', 'Bex2', 'Rph3al', 'Scg5', 'Pcsk1', 'Isl1', 'Maged1', 'Fabp5', 'Celf3', 'Pcsk1n', 'Fam183b', 'Prnp', 'Tac1', 'Gpx3', 'Cplx2', 'Nkx2-2', 'Olfm1', 'Vim', 'Rimbp2', 'Anxa6', 'Scg3', 'Ngfrap1', 'Insm1', 'Gng4', 'Pax6', 'Cnot6l', 'Cacna2d1', 'Tox3', 'Slc39a2', 'Riiad1']
    marker_genes['Tuft'] = ['Alox5ap', 'Lrmp', 'Hck', 'Avil', 'Rgs13', 'Ltc4s', 'Trpm5', 'Dclk1', 'Spib', 'Fyb', 'Ptpn6', 'Matk', 'Snrnp25', 'Sh2d7', 'Ly6g6f', 'Kctd12', '1810046K07Rik', 'Hpgds', 'Tuba1a', 'Pik3r5', 'Vav1', 'Tspan6', 'Skap2', 'Pygl', 'Ccdc109b', 'Ccdc28b', 'Plcg2', 'Ly6g6d', 'Alox5', 'Pou2f3', 'Gng13', 'Bmx', 'Ptpn18', 'Nebl', 'Limd2', 'Pea15a', 'Tmem176a', 'Smpx', 'Itpr2', 'Il13ra1', 'Siglecf', 'Ffar3', 'Rac2', 'Hmx2', 'Bpgm', 'Inpp5j', 'Ptgs1', 'Aldh2', 'Pik3cg', 'Cd24a', 'Ethe1', 'Inpp5d', 'Krt23', 'Gprc5c', 'Reep5', 'Csk', 'Bcl2l14', 'Tmem141', 'Coprs', 'Tmem176b', '1110007C09Rik', 'Ildr1', 'Galk1', 'Zfp428', 'Rgs2', 'Inpp5b', 'Gnai2', 'Pla2g4a', 'Acot7', 'Rbm38', 'Gga2', 'Myo1b', 'Adh1', 'Bub3', 'Sec14l1', 'Asah1', 'Ppp3ca', 'Agt', 'Gimap1', 'Krt18', 'Pim3', '2210016L21Rik', 'Tmem9', 'Lima1', 'Fam221a', 'Nt5c3', 'Atp2a3', 'Mlip', 'Vdac3', 'Ccdc23', 'Tmem45b', 'Cd47', 'Lect2', 'Pla2g16', 'Mocs2', 'Arpc5', 'Ndufaf3']
    cell_annotation = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r0.5')
    cell_annotation_norm = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r0.5', normalize='reference')
    sb.heatmap(cell_annotation_norm, cbar=False, annot=True)

    # Define a nice colour map for gene expression
    colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
    colors3 = plt.cm.Greys_r(np.linspace(0.7, 0.8, 20))
    colorsComb = np.vstack([colors3, colors2])
    mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
    # Defa24 #Tff3
    sc.pl.umap(adata, color='Defa24', use_raw=False, color_map=mymap)
    sc.pl.umap(adata, color='Tff3', use_raw=False, color_map=mymap)

    # Check expression of enterocyte markers
    # Collate all enterocyte markers and get the gene IDs in the data set
    ids_entprox = np.in1d(adata.var_names, marker_genes['Enterocyte (Proximal)'])
    ids_entdist = np.in1d(adata.var_names, marker_genes['Enterocyte (Distal)'])
    ids_ent = np.logical_or(ids_entprox, ids_entdist)

    # Calculate the mean expression of enterocyte markers
    adata.obs['Enterocyte_marker_expr'] = adata.X[:, ids_ent].mean(1)

    # Plot enterocyte expression
    sc.pl.violin(adata, 'Enterocyte_marker_expr', groupby='louvain_r0.5')
    sc.pl.umap(adata, color='Enterocyte_marker_expr', color_map=mymap)

    #Early enterocyte marker - Arg2
    sc.pl.umap(adata, color='Arg2', use_raw=False, color_map=mymap)

    sc.pl.violin(adata, groupby='louvain_r0.5', keys='Arg2', use_raw=False)

    sc.pl.diffmap(adata, components=['6,9'], color='Arg2', use_raw=False, color_map=mymap)
    sc.pl.diffmap(adata, components=['6,9'], color='louvain_r0.5')

    sc.pl.violin(adata, 'mt_frac', groupby='louvain_r0.5')
    sc.pl.violin(adata, 'log_counts', groupby='louvain_r0.5')

    # Check individual stem markers
    stem_genes = adata.var_names[np.in1d(adata.var_names, marker_genes['Stem'])]
    sc.pl.umap(adata, color=stem_genes[:3], title=stem_genes[:3], color_map=mymap)
    sc.pl.umap(adata, color=stem_genes[3:], title=stem_genes[3:], color_map=mymap)

    # Check stem marker expression
    adata.obs['Stem_marker_expr'] = adata[:, stem_genes].X.mean(1)

    sc.pl.violin(adata, 'Stem_marker_expr', groupby='louvain_r0.5')
    sc.pl.umap(adata, color='Stem_marker_expr', color_map=mymap)

    # Categories to rename
    adata.obs['louvain_r0.5'].cat.categories
    adata.rename_categories('louvain_r0.5', ['TA', 'EP (early)', 'Stem', 'Goblet', 'EP (stress)', 'Enterocyte', 'Paneth', 'Enteroendocrine', 'Tuft'])

    adata.obs['louvain_r0.5'].value_counts()

    sc.pl.umap(adata, color='louvain_r0.5', size=15, legend_loc='on data')




def main():
    sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
    cc_genes_file = '~/data/Macosko_cell_cycle_genes.txt'
    # cc_genes_file = 'E:\\sc_data\\Macosko_cell_cycle_genes.txt'

    # 获得数据
    adata = get_raw_data()
    adata = pre_qc(adata)
    #
    # # qc 画图
    # ## 查看在样本中的分布
    # sc_qc_violin_plot("n_counts", group='sample', **{"size": 1, "log": True, "show": False})
    # sc_qc_violin_plot("mt_frac", group='sample', **{"size": 1, "show": False})
    #
    # ## 查看数据之间的关系
    # sc_qc_scatter_plot(adata, 'n_counts', 'n_genes', color='mt_frac', save=1)
    # sc_qc_scatter_plot(adata[adata.obs['n_counts'] < 10000], 'n_counts', 'n_genes', color='mt_frac', save=2)
    #
    # ## 质控，确定删选的阈值
    # sc_qc_dist_plot(adata.obs['n_counts'], "n_counts", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['n_counts'][adata.obs['n_counts'] < 4000], "n_counts_less_4000", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['n_counts'][adata.obs['n_counts'] > 10000], "n_counts_more_10000", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['n_genes'], "n_genes", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['n_genes'][adata.obs['n_genes'] < 1000], "n_genes_less_1000", kde=False, bins=60)
    # sc_qc_dist_plot(adata.obs['mt_frac'], "mt_frac", kde=False, bins=60)
    # # 以上得到n_counts范围1500~40000；n_genes范围为>700;mt_frac范围为<0.2
    #
    # filter cells
    adata = sc_qc_filter(adata, min_counts=1500, max_counts=40000, min_genes=700, mt_frac=0.2)
    #
    # Filter genes:
    print('Total number of genes: {:d}'.format(adata.n_vars))
    # Min 20 cells - filters out 0 count genes
    sc.pp.filter_genes(adata, min_cells=20)
    print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
    #
    #
    # # select top genes
    # sc_top_gene_plot(adata, sample_strings)


    # normalization

    # # TODO 解释以下代码意义 --DONE
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)  # 百万细胞
    sc.pp.log1p(adata_pp)   # 以自然数为底，log转换
    sc.pp.pca(adata_pp, n_comps=15)  # PCA分析
    sc.pp.neighbors(adata_pp)
    sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)  # 聚类分析
    # # Normalization 得到size_factors
    # size_factors = get_size_factors(adata_pp)
    # del adata_pp
    # size_factors.index = adata.obs.index
    # adata.obs["size_factors"] = size_factors
    #
    # sc_qc_scatter_plot(adata, 'size_factors', 'n_counts',save="size_factors_n_counts")
    # sc_qc_scatter_plot(adata, 'size_factors', 'n_genes',save="size_factors_n_genes")
    # sc_qc_dist_plot(adata.obs['size_factors'],'size_factors', bins=50, kde=False)
    #
    # # Keep the count data in a counts layer
    # adata.layers["counts"] = adata.X.copy()
    # # Normalize adata
    # adata.X /= adata.obs['size_factors'].values[:, None]
    # sc.pp.log1p(adata)
    #
    # # Batch Correction
    # # sc.pp.combat(adata, key='sample')
    #
    # # Highly Variable Genes
    # sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
    # print('\n', 'Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
    # sc.pl.highly_variable_genes(adata,show=False, save="highly_variable_genes.png")
    #
    # to_visualization(adata)
    # cell_cycle_scoring(adata,cc_genes_file)
    # to_clustering(adata)


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    main()






