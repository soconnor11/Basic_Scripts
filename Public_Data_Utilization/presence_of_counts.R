


write.csv(rowSums(data@assays$RNA@counts>=3)/ncol(data), file.path(resdir, paste0(tag, "_presence_gte3.csv")))
