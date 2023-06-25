library(ggrepel)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(UCell)
library(SingleR)
library(limma)
library(stringr)
library(jsonlite) # 用于将结果转为json
library(org.Hs.eg.db)
library(patchwork)
library(presto)
library(scRNAtoolVis)
library(corrplot)

### 方法内判断输入样本的格式
readData <- function(dirPath, fileList, project) {
  filesCount = length(fileList)
  if (filesCount == 0) {
    print("样本路径读取为空")
    return(NULL)
  }else if (filesCount == 1) {
    obj.counts <- NULL
    if (endsWith(fileList, ".csv")) {
      obj.counts = read.csv(fileList, header = T, row.names = 1)
    }else if (endsWith(fileList, ".txt")) {
      obj.counts = read.table(fileList, header = T)
    }
    print(fileList)
    obj = CreateSeuratObject(obj.counts, project = project, min.cells = 3, min.features = 200)
    obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^MT-")
    return(obj)
  }else if (filesCount == 3) {
    print(dirPath)
    obj.counts = Read10X(dirPath)
    obj = CreateSeuratObject(obj.counts, project = project, min.cells = 3, min.features = 200)
    obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^MT-")
    return(obj)
  }
  return(NULL)
}

### 判断输入样本的格式
dir <- "/data/sca/input/202306/22/uid1_20230622014435_ef0ab46794de477bb32a670bfbbef4fb"
fileList <- list.files(dir, include.dirs = F, full.names = TRUE, recursive = F)
pbmc_expr = readData(dir, fileList, "expr")

jsonlite::toJSON(data.frame(gene = dim(pbmc_expr)[1], cell = dim(pbmc_expr)[2]), pretty = F)

dir <- "/data/sca/control_group_data/30"
fileList <- list.files(dir, include.dirs = F, full.names = TRUE, recursive = F)
pbmc_ctrl = readData(dir, fileList, "ctrl")

pbmc = merge(pbmc_expr, pbmc_ctrl, add.cell.ids = c("expr", "ctrl"))

## 添加分组信息
pbmc$sample = stringr::str_split_fixed(colnames(pbmc), "_", n = 2)[, 1]


### 质量控制 QC

nFeature_RNA_value <- round(as.matrix(quantile(pbmc$nFeature_RNA, 96 / 100))[1], 2)
nCount_RNA_value <- round(as.matrix(quantile(pbmc$nCount_RNA, 96 / 100))[1], 2)
percent_mt_value <- round(as.matrix(quantile(pbmc$percent.mt, 90 / 100))[1], 2)

p1 <- VlnPlot(pbmc, features = "percent.mt") &
  geom_hline(linetype = 'dotdash', col = 'red', yintercept = percent_mt_value, size = 1) &
  NoLegend() &
  annotate(geom = "label", x = 2, y = percent_mt_value, label = percent_mt_value)
p2 <- VlnPlot(pbmc, features = "nCount_RNA") &
  geom_hline(linetype = 'dotdash', col = 'red', yintercept = nCount_RNA_value, size = 1) &
  NoLegend() &
  annotate(geom = "label", x = 2, y = nCount_RNA_value, label = nCount_RNA_value)
p3 <- VlnPlot(pbmc, features = "nFeature_RNA") &
  geom_hline(linetype = 'dotdash', col = 'red', yintercept = nFeature_RNA_value, size = 1) &
  NoLegend() &
  annotate(geom = "label", x = 2, y = nFeature_RNA_value, label = nFeature_RNA_value)

pQC <- wrap_plots(p1, p2, p3, ncol = 3)

# pQC<-VlnPlot(pbmc, features = c("percent.mt","nCount_RNA","nFeature_RNA"), ncol = 3,pt.size =0.1) ####查看数据原始分布情况

ggsave(
  filename = "/data/sca/user_data/18/output/plot_qc.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 2000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 300,              # 分辨率DPI
  plot = pQC,
  limitsize = FALSE
)


#### 通过一定的指标来进行过滤的选择
# aging <- subset(pbmc,nCount_RNA >= 800 & nCount_RNA<tail(hist(pbmc$nCount_RNA)$mids,1) & percent.mt <=30 & nFeature_RNA <tail(hist(pbmc$nFeature_RNA)$mids,1) & nFeature_RNA>500)
aging <- subset(pbmc, nCount_RNA >= 800 &
  nCount_RNA < nCount_RNA_value &
  percent.mt <= percent_mt_value &
  nFeature_RNA < nFeature_RNA_value &
  nFeature_RNA > 500)

### 归一化后pca降维，寻找合适的维度拐点
# 归一化
aging <- NormalizeData(aging, normalization.method = "LogNormalize", scale.factor = 10000) ####默认参数
# 寻找高变异基因->scale归一化->跑pca降维(主成分分析)
aging <- FindVariableFeatures(aging, selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA()  ###nfeatures一般选2000-5000，对结果影响较大，需要手动选择
# 从拐点图选择合适的维度值
pca_num <- ElbowPlot(aging, ndims = 40)

ggsave(
  filename = "/data/sca/user_data/18/output/pca_dim_num.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 2000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 300,              # 分辨率DPI
  plot = pca_num,
  limitsize = FALSE
)


# 去批次
library(harmony)
aging <- RunHarmony(aging, group.by.vars = "sample")
#降维聚类
aging <- FindNeighbors(aging, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.2)
aging <- RunUMAP(aging, reduction = "harmony", dims = 1:20, label = T) %>% RunTSNE(reduction = "harmony", dims = 1:20, label = T)

## >>>> 出图
pumap <- DimPlot(aging, reduction = "umap", group.by = c("sample"), label = T)

ggsave(
  filename = "/data/sca/user_data/18/output/plot_umap.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 2600,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 300,              # 分辨率DPI
  plot = pumap,
  limitsize = FALSE
)


#aging.markers=FindAllMarkers(aging,only.pos = T,assay = "RNA",logfc.threshold = 0.25)
aging.markers <- subset(wilcoxauc(aging, "seurat_clusters"), logFC >= 0.15 & pct_in >= 0.1 & pct_out >= 0.1) %>% dplyr::rename(gene = feature, cluster = group, avg_log2FC = logFC)


# rda=>load rds=>read
#load("/data/sca/ref/reference.rds")
#load("/data/sca/ref/20230530-test_Lu.rda")
#load("/data/sca/ref/scRNA_zhushilabel_SE.ref2.rds")
young1 <- readRDS("/data/sca/ref/scRNA_zhushi.rds")

test_label <- t(FetchData(young1, vars = c("ident")))

test_data <- as.data.frame(GetAssayData(young1, slot = "data"))
print("test")

test_ref_list <- list(count = test_data, label = test_label)
print("test")

aging_for_SingleR <- GetAssayData(aging, slot = "data") ##获取标准化矩阵
aging.hesc <- SingleR(test = aging_for_SingleR, ref = test_ref_list$count, labels = test_ref_list$label)
print("test")

aging@meta.data$labels <- aging.hesc$labels

## 将注释的label加到ident中
Idents(aging) <- aging$labels
## 定义细胞类型
aging$celltype <- aging@active.ident

plot_celltytpe <- DimPlot(aging, group.by = c("seurat_clusters", "labels"), reduction = "umap", label = T)

ggsave(
  filename = "/data/sca/user_data/18/output/annotation_result.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = plot_celltytpe,
  limitsize = FALSE
)

### 测试 Marker 表达  出图
genes_to_check <- c("DAZL", "DDX4", "MAGEA4", "UTF1", "FCGR3A", "KIT", "DMRT1", "DMRTB1", "STRA8", "SYCP3", "SPO11", "MLH3", "ZPBP", "ID4", "PIWIL4", "UCHL1", "TNP1", "TNP2", "PRM2", "SOX9", "WT1", "AMH", "PRND", "FATE1", "VWF", "PECAM1", "CDH5", "DLK1", "IGF1", "CYP11A1", "STAR", "NOTCH3", "ACTA2", "MYH11", "CYP26B1", "WFDC1", "CD14", "CD163", "C1QA", "C1QC", "CD8A", "CD8B", "PTPRC")
p_all_markers <- DotPlot(aging, features = genes_to_check, assay = 'RNA', group.by = 'celltype') + coord_flip()

dot_data <- p_all_markers$data

colnames(dot_data) <- c("AverageExpression_unscaled", "Precent Expressed", "Features", "celltype", "Average Expression")

####用ggplot画图####
p_marker_dotplot = ggplot(dot_data, aes(celltype, Features, size = `Precent Expressed`)) +
  geom_point(shape = 21, aes(fill = `Average Expression`), position = position_dodge(0)) +
  theme_minimal() +
  xlab(NULL) +
  ylab(NULL) +
  scale_size_continuous(range = c(1, 10)) +
  theme_bw() +
  scale_fill_gradient(low = "grey", high = "#E54924") +
  theme(legend.position = "right", legend.box = "vertical", #图例位置
        legend.margin = margin(t = 0, unit = 'cm'),
        legend.spacing = unit(0, "in"),
        axis.text.x = element_text(color = "black", size = 16, angle = 45,
                                   hjust = 1), #x轴
        axis.text.y = element_text(color = "black", size = 12), #y轴
        legend.text = element_text(size = 12, color = "black"), #图例
        legend.title = element_text(size = 12, color = "black"), #图例
        axis.title.y = element_text(vjust = 1,
                                    size = 16)
  ) +
  labs(x = " ", y = "Features")

ggsave(
  filename = "E:/singlecelltest/test_project/GSE112013_Combined_UMI_table/test_data_output/p_marker_dotplot.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = p_marker_dotplot,
  limitsize = FALSE
)


## 细胞比例计算
## 样本间的比例变化
celltype_ratio <- prop.table(table(Idents(aging), aging$sample), margin = 2)
celltype_ratio <- as.data.frame(celltype_ratio)
colnames(celltype_ratio)[1] <- "celltype"
colnames(celltype_ratio)[2] <- "sample"
colnames(celltype_ratio)[3] <- "freq"
## celltype_ratio是个数据框，可以提取细胞比例信息

p_ratio <- ggplot(celltype_ratio) +
  geom_bar(aes(x = freq, y = sample, fill = celltype), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = 'Ratio', y = 'Sample') +
  coord_flip() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

jsonlite::toJSON(celltype_ratio, pretty = F)

ggsave(
  filename = "/data/sca/user_data/18/output/p_ratio.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = p_ratio,
  limitsize = FALSE
)


## 分开生殖细胞与体细胞
## 生殖细胞
germ_cell <- c("Round S'tids", "Sperm", "Elongated S'tids", "Late primary S'gonia", "SSCs", "Early primary S'gonia", "Differentiating S'gonia")
germ <- subset(aging, idents = germ_cell)
germ$celltype <- germ@active.ident

## 体细胞
somatic <- subset(aging, idents = germ_cell, invert = TRUE)
somatic$celltype <- somatic@active.ident
## 生殖细胞UMAP


## 定义运行的函数
UmapCellratioFun <- function(cellOBj, celltype_prefix) {
  ### 归一化后pca降维，寻找合适的维度拐点
  # 归一化
  cellOBj <- NormalizeData(cellOBj, normalization.method = "LogNormalize", scale.factor = 10000) ####默认参数
  # 寻找高变异基因->scale归一化->跑pca降维(主成分分析)
  cellOBj <- FindVariableFeatures(cellOBj, selection.method = "vst", nfeatures = 1000) %>%
    ScaleData() %>%
    RunPCA()  ###nfeatures一般选2000-5000，对结果影响较大，需要手动选择
  # 从拐点图选择合适的维度值
  pca_num <- ElbowPlot(cellOBj, ndims = 40)

  # 去批次
  library(harmony)
  cellOBj <- RunHarmony(cellOBj, group.by.vars = "sample")
  #降维聚类
  cellOBj <- FindNeighbors(cellOBj, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.2)
  cellOBj <- RunUMAP(cellOBj, reduction = "harmony", dims = 1:20, label = T) %>% RunTSNE(reduction = "harmony", dims = 1:20, label = T)

  ## 用整体注释的细胞类型重新定义生殖细胞
  Idents(cellOBj) <- cellOBj$celltype

  ## >>>> 出图
  pumap <- DimPlot(cellOBj, reduction = "umap")

  ggsave(
    filename = paste0("/data/sca/user_data/18/output/", celltype_prefix, "UMAP.png", seq = ""), # 保存的文件名称。通过后缀来决定生成什么格式的图片
    width = 2600,             # 宽
    height = 2000,            # 高
    units = "px",          # 单位
    dpi = 300,              # 分辨率DPI
    plot = pumap,
    limitsize = FALSE
  )

  ## 细胞比例计算
  ## 样本间的比例变化
  celltype_ratio <- prop.table(table(Idents(cellOBj), cellOBj$sample), margin = 2)
  celltype_ratio <- as.data.frame(celltype_ratio)
  colnames(celltype_ratio)[1] <- "celltype"
  colnames(celltype_ratio)[2] <- "sample"
  colnames(celltype_ratio)[3] <- "freq"
  ## celltype_ratio是个数据框，可以提取细胞比例信息

  p_ratio <- ggplot(celltype_ratio) +
    geom_bar(aes(x = freq, y = sample, fill = celltype), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') +
    theme_classic() +
    labs(x = 'Ratio', y = 'Sample') +
    coord_flip() +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

  jsonlite::toJSON(celltype_ratio, pretty = F)

  ggsave(
    filename = paste0("/data/sca/user_data/18/output/", celltype_prefix, "CellRatio.png", seq = ""), # 保存的文件名称。通过后缀来决定生成什么格式的图片
    width = 4000,             # 宽
    height = 2000,            # 高
    units = "px",          # 单位
    dpi = 250,              # 分辨率DPI
    plot = p_ratio,
    limitsize = FALSE
  )

  ## 差异分析

  # 基因差异分析
  DEG_all <- rbind.data.frame()

  for (item in rev(unique(cellOBj$celltype))) {
    tmp <- subset(cellOBj, idents = item)
    errCheck = tryCatch({
      tmp.markers <- FindMarkers(tmp, group.by = "sample", ident.1 = "expr", ident.2 = "ctrl", min.pct = 0.1, logfc.threshold = 0.25)
      tmp.markers$gene <- rownames(tmp.markers)
      0
    }, error = function(e) {
      # print(e)
      2
    })
    if (errCheck == 2) {
      # print("遇到错误，结束循环")
      next
    }
    tmp.markers$cluster <- item
    DEG_all <- rbind.data.frame(DEG_all, tmp.markers)
  }

  write.csv(DEG_all, file = paste0("/data/sca/user_data/18/output/", celltype_prefix, "different_expression_gene.csv", seq = ""))

  p_deg_volcano <- jjVolcano(diffData = DEG_all, tile.col = corrplot::COL2('RdBu', length(table(DEG_all$cluster))))

  ggsave(
    filename = paste0("/data/sca/user_data/18/output/", celltype_prefix, "p_deg_volcano.png", seq = ""),
    width = 4000,             # 宽
    height = 2000,            # 高
    units = "px",          # 单位
    dpi = 250,              # 分辨率DPI
    plot = p_deg_volcano,
    limitsize = FALSE
  )

  ####差异基因数统计
  count = subset(DEG_all, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
  count = data.frame(table(count$cluster, count$avg_log2FC > 0.25))
  colnames(count) = c("cluster", "log2FC", "Freq")
  count$Freq2 = ifelse(count$log2FC == "TRUE", count$Freq, 0 - count$Freq)
  count$fill = ifelse(count$log2FC == "TRUE", "Up", "Down")
  count$fill = factor(count$fill, levels = c("Up", "Down"))
  deg_bar_plot = ggplot(count, aes(x = cluster, y = Freq2, fill = fill)) +
    geom_bar(stat = 'identity', position = 'stack') +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 13, angle = 0, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 13),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15)) +
    guides(fill = guide_legend(title = NULL)) +
    geom_text(label = count$Freq, nudge_x = 0, nudge_y = 1) +
    xlab("Clusters") +  #x轴标签
    ylab("DEG counts") +  #y轴标签
    labs(title = "Differention Expression Gene Counts")  #设置标题

  ggsave(
    filename = paste0("/data/sca/user_data/18/output/", celltype_prefix, "deg_bar_plot.png", seq = ""),
    width = 4000,             # 宽
    height = 2000,            # 高
    units = "px",          # 单位
    dpi = 250,              # 分辨率DPI
    plot = deg_bar_plot,
    limitsize = FALSE)


  library(clusterProfiler)

  #### 针对差异基因进行通路富集分析，区分上调基因和下调基因
  go.enrich = function(gene) {
    eg = bitr(gene, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb = "org.Hs.eg.db")
    ego <- enrichGO(gene = eg[, 2], OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.4,
                    qvalueCutoff = 0.2, readable = T)
    if (is.null(ego) ||
      is.null(ego@result) ||
      length(rownames(ego@result)) == 0) {
      return(NULL);
    }
    go = data.frame(ego@result)
    go$GeneRatio2 <- sapply(go$GeneRatio, function(x) eval(parse(text = x)))
    return(go)
  }

  ### 细胞类型的富集分析,针对基因的
  gene.enrich = function(data) {
    go <- rbind.data.frame()
    for (i in unique(data$cluster)) {
      test = subset(data, cluster == i)
      result = go.enrich(test$gene)
      if (is.null(result)) {
        next;
      }
      result$cluster = i
      go = rbind(go, result)
    }
    return(go)
  }


  ### 上调基因
  up = subset(DEG_all, p_val_adj < 0.05 & avg_log2FC > 0.25 )
  up.go = gene.enrich(up)
  #####筛选通路
  # 定义要排除的关键词
  keywords <- c("cytoplasmic translation","aerobic respiration","ribosome,rRNA processing","rRNA metabolic process","translational initiation","RNA capping","ncRNA processing","regulation of translation","RNA splicing","mRNA splicing","protein targeting","translational elongation","maturation of SSU-rRNA","respiratory burst","entry into host","fructose,maturation of LSU-rRNA","ribosome localization","cytochrome complex assembly","'de novo' protein folding","lung alveolus development","ossification","artery morphogenesis","lung development","ovulation","regulation of body fluid levels","liver morphogenesis","lung morphogenesis","neural tube development","neural tube closure","endoderm formation","exploration behavior","tube closure","lung cell differentiation","lung epithelial cell differentiation","digestive tract development","primary neural tube formation","hair follicle morphogenesis","RNA modification","vesicle cargo loading","virion assembly","ovulation cycle","response to food","locomotory behavior","Golgi vesicle budding","centriole replication","bleb assembly","biological phase","elastic fiber assembly","lactation","luteinization","visual system development","endoderm development","telencephalon development","reproductive behavior","ruffle assembly","centriole assembly","associative learning","adhesion of symbiont to host cell","visual learning","intracellular lipid transport","protein destabilization","mitochondrial translation","cellular response to starvation","response to caffeine","radial glial cell differentiation","neuron migration","centriole-centriole cohesion","motor behavior","protein homotrimerization","left/right axis specification","mitochondrion disassembly","axis specification","response to dsRNA","trophoblast cell migration","glial cell migration","oocyte development","embryo implantation","mRNA cleavage","mRNA stabilization","mRNA catabolic process","hindbrain development","oocyte differentiation","sialylation","heart contraction","RNA catabolic process","heart looping","protein sumoylation","cardiocyte differentiation","eye morphogenesis","histone H2A monoubiquitination","heart process","nucleic acid transport","RNA transport","mitochondrial gene expression","maturation of 5.8S rRNA from tricistronic rRNA transcript","translation","translational initiation","RNA binding","RNA splicing","protein binding","transcription","RNA polymerase II transcriptional preinitiation complex assembly")
  # 使用grepl()函数进行筛选
  up.go <- up.go[!grepl(paste(keywords, collapse = "|"), up.go$Description), ]
  up.go <- subset(up.go, Count >= 2)
  write.csv(up.go, file = paste0("/data/sca/user_data/18/output/", celltype_prefix, "upGo.csv", seq = ""))

  upGoTopN <- up.go %>%
    group_by(cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = 15) %>%
    arrange(desc(Count))

  ### 下调基因
  down = subset(DEG_all, p_val_adj < 0.05 & avg_log2FC < -0.25)
  down.go = gene.enrich(down)
  #####筛选通路
  # 定义要排除的关键词
  keywords <- c("cytoplasmic translation","aerobic respiration","ribosome,rRNA processing","rRNA metabolic process","translational initiation","RNA capping","ncRNA processing","regulation of translation","RNA splicing","mRNA splicing","protein targeting","translational elongation","maturation of SSU-rRNA","respiratory burst","entry into host","fructose,maturation of LSU-rRNA","ribosome localization","cytochrome complex assembly","'de novo' protein folding","lung alveolus development","ossification","artery morphogenesis","lung development","ovulation","regulation of body fluid levels","liver morphogenesis","lung morphogenesis","neural tube development","neural tube closure","endoderm formation","exploration behavior","tube closure","lung cell differentiation","lung epithelial cell differentiation","digestive tract development","primary neural tube formation","hair follicle morphogenesis","RNA modification","vesicle cargo loading","virion assembly","ovulation cycle","response to food","locomotory behavior","Golgi vesicle budding","centriole replication","bleb assembly","biological phase","elastic fiber assembly","lactation","luteinization","visual system development","endoderm development","telencephalon development","reproductive behavior","ruffle assembly","centriole assembly","associative learning","adhesion of symbiont to host cell","visual learning","intracellular lipid transport","protein destabilization","mitochondrial translation","cellular response to starvation","response to caffeine","radial glial cell differentiation","neuron migration","centriole-centriole cohesion","motor behavior","protein homotrimerization","left/right axis specification","mitochondrion disassembly","axis specification","response to dsRNA","trophoblast cell migration","glial cell migration","oocyte development","embryo implantation","mRNA cleavage","mRNA stabilization","mRNA catabolic process","hindbrain development","oocyte differentiation","sialylation","heart contraction","RNA catabolic process","heart looping","protein sumoylation","cardiocyte differentiation","eye morphogenesis","histone H2A monoubiquitination","heart process","nucleic acid transport","RNA transport","mitochondrial gene expression","maturation of 5.8S rRNA from tricistronic rRNA transcript","translation","translational initiation","RNA binding","RNA splicing","protein binding","transcription","RNA polymerase II transcriptional preinitiation complex assembly")
  # 使用grepl()函数进行筛选
  down.go <- down.go[!grepl(paste(keywords, collapse = "|"), down.go$Description), ]
  down.go <- subset(down.go, Count >= 2)
  write.csv(down.go, file = paste0("/data/sca/user_data/18/output/", celltype_prefix, "downGo.csv", seq = ""))

  downGoTopN <- down.go %>%
    group_by(cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = 15) %>%
    arrange(desc(Count))


  #点图#
  up_go_point_plot <- ggplot(upGoTopN, aes(x = cluster, y = reorder(Description, -pvalue), size = Count, color = -log10(pvalue))) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 13, angle = 0, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 13),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15)) +
    scale_color_gradient(low = "lightgrey", high = "red") +
    xlab("Clusters") +  #x轴标签
    ylab("Pathway") +  #y轴标签
    labs(title = "Up Regulate GO Terms Enrichment") +  #设置标题
    facet_wrap(~ONTOLOGY, ncol = 3)

  ggsave(
    filename = paste0("/data/sca/user_data/18/output/", celltype_prefix, "up_go_point_plot.png", seq = ""),
    width = 4000,             # 宽
    height = 2000,            # 高
    units = "px",          # 单位
    dpi = 250,              # 分辨率DPI
    plot = up_go_point_plot,
    limitsize = FALSE
  )


  # 下调基因
  #纵向点图#
  down_go_point_plot <- ggplot(downGoTopN, aes(x = cluster, y = reorder(Description, -pvalue), size = Count, color = -log10(pvalue))) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 13, angle = 0, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 13),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    xlab("Clusters") +  #x轴标签
    ylab("Pathway") +  #y轴标签
    labs(title = "Down Regulate GO Terms Enrichment") +  #设置标题
    facet_wrap(~ONTOLOGY, ncol = 3)

  ggsave(
    filename = paste0("/data/sca/user_data/18/output/", celltype_prefix, "down_go_point_plot.png", seq = ""),
    width = 4000,             # 宽
    height = 2000,            # 高
    units = "px",          # 单位
    dpi = 250,              # 分辨率DPI
    plot = down_go_point_plot,
    limitsize = FALSE
  )

}

## 调用函数
## 运行生殖细胞
UmapCellratioFun(germ, "germ")

## 运行体细胞
UmapCellratioFun(somatic, "somatic")


