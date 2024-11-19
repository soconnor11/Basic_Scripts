
SaveH5Seurat(data, filename = "hBMMSC_ATCC.h5Seurat")
Convert("hBMMSC_ATCC.h5Seurat", dest = "h5ad")


adata_vc = scv.read_loom('hBMMSC_CKDL210009543-1a-SI_GA_G1_H57KLDSX2/velocyto/hBMMSC_CKDL210009543-1a-SI_GA_G1_H57KLDSX2.loom')
adata_Seurat = sc.read_h5ad('redo_analysis/hBMMSC_ATCC.h5ad')
adata = scv.utils.merge(adata_vc, adata_Seurat)
adata.write_h5ad("hBMMSC_ATCC_with_velocyto_info.h5ad")
