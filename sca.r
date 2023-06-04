
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

# 读取样本数据 10x方式 3个文件的方式
creatobj=function(file,project){
  obj.counts=Read10X(file)
  obj=CreateSeuratObject(obj.counts,project = project,min.cells=3,min.features=200)
  obj[["percent.mt"]]=PercentageFeatureSet(obj, pattern = "^MT-")
  return(obj)}

# 读取样本数据 csv方式  1个文件的方式
creatobj2=function(file,project){
  obj.counts=read.table(file,header = T)
  obj=CreateSeuratObject(obj.counts,project = project,min.cells=3,min.features=200)
  obj[["percent.mt"]]=PercentageFeatureSet(obj, pattern = "^MT-")
  return(obj)}

### 判断输入样本的格式
dir <- "/data/sca/input/202306/02/uid1_20230602220441_69c150b9d0464f9c90865f9543cd8014"
flist <- list.files(dir, include.dirs = F,full.names = TRUE,recursive = F)

### create S4格式
if(length(flist)==3){
  print("加载实验组3个文件")
  pbmc_expr=creatobj(dir,"expr")
}else{
  print("加载实验组1个文件")
  pbmc_expr=creatobj2(paste0(dir,"/*",collapse = NULL),"expr") ###这里project可以换成疾病的类型
}

jsonlite::toJSON(data.frame(gene = dim(pbmc_expr)[1], cell = dim(pbmc_expr)[2]), pretty = F)

### 判断输入样本的格式
dir <- "/data/sca/control_group_data/30"
flist <- list.files(dir, include.dirs = F,full.names = TRUE,recursive = F)

### create S4格式
if(length(flist)==3){
  print("加载对照组3个文件")
  pbmc_ctrl=creatobj(dir,"ctrl")
}else{
  print("加载对照组1个文件")
  pbmc_expr=creatobj2(paste0(dir,"/*",collapse = NULL),"ctrl") ###这里project可以换成疾病的类型
}

pbmc=merge(pbmc_expr,pbmc_ctrl,add.cell.ids = c("expr","ctrl"))

## 添加分组信息
pbmc$sample=stringr::str_split_fixed(colnames(pbmc),"_",n=2)[,1]


### 质量控制 QC

nFeature_RNA_value <- round(as.matrix(quantile(pbmc$nFeature_RNA,96/100))[1],2)
nCount_RNA_value <- round(as.matrix(quantile(pbmc$nCount_RNA,96/100))[1],2)
percent_mt_value <- round(as.matrix(quantile(pbmc$percent.mt,90/100))[1],2)

p1 <- VlnPlot(pbmc, features = "percent.mt") & geom_hline(linetype='dotdash', col = 'red', yintercept = percent_mt_value, size=1) & NoLegend() & annotate(geom = "label", x=2, y= percent_mt_value, label=percent_mt_value)
p2 <- VlnPlot(pbmc, features = "nCount_RNA") & geom_hline(linetype='dotdash', col = 'red', yintercept = nCount_RNA_value, size=1) & NoLegend() & annotate(geom = "label", x=2, y= nCount_RNA_value, label=nCount_RNA_value)
p3 <- VlnPlot(pbmc, features = "nFeature_RNA") & geom_hline(linetype='dotdash', col = 'red', yintercept = nFeature_RNA_value, size=1) & NoLegend() & annotate(geom = "label", x=2, y= nFeature_RNA_value, label=nFeature_RNA_value)

pQC <- wrap_plots(p1, p2, p3, ncol = 3)

# pQC<-VlnPlot(pbmc, features = c("percent.mt","nCount_RNA","nFeature_RNA"), ncol = 3,pt.size =0.1) ####查看数据原始分布情况

ggsave(
  filename = "/data/sca/user_data/4/output/plot_qc.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 2000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 300,              # 分辨率DPI
  plot = pQC,
  limitsize = FALSE
)


#### 通过一定的指标来进行过滤的选择
# aging <- subset(pbmc,nCount_RNA >= 800 & nCount_RNA<tail(hist(pbmc$nCount_RNA)$mids,1) & percent.mt <=30 & nFeature_RNA <tail(hist(pbmc$nFeature_RNA)$mids,1) & nFeature_RNA>500)
aging <- subset(pbmc,nCount_RNA >= 800 & nCount_RNA< nCount_RNA_value & percent.mt <= percent_mt_value & nFeature_RNA < nFeature_RNA_value & nFeature_RNA>500)

### 归一化后pca降维，寻找合适的维度拐点
# 归一化
aging <- NormalizeData(aging, normalization.method = "LogNormalize", scale.factor = 10000)####默认参数
# 寻找高变异基因->scale归一化->跑pca降维(主成分分析)
aging <- FindVariableFeatures(aging, selection.method = "vst", nfeatures = 2000)%>%ScaleData()%>%RunPCA()  ###nfeatures一般选2000-5000，对结果影响较大，需要手动选择
# 从拐点图选择合适的维度值
pca_num<-ElbowPlot(aging, ndims = 40)

ggsave(
  filename = "/data/sca/user_data/4/output/pca_dim_num.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
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
aging <- RunUMAP(aging, reduction = "harmony", dims = 1:20,label = T) %>% RunTSNE(reduction = "harmony",dims = 1:20,label = T)

## >>>> 出图
pumap<-DimPlot(aging, reduction = "umap",group.by=c("sample"),label = T)

ggsave(
  filename = "/data/sca/user_data/4/output/plot_umap.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 2600,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 300,              # 分辨率DPI
  plot = pumap,
  limitsize = FALSE
)


#aging.markers=FindAllMarkers(aging,only.pos = T,assay = "RNA",logfc.threshold = 0.25)
aging.markers <- subset(wilcoxauc(aging, "seurat_clusters"), logFC>=0.15&pct_in >=0.1&pct_out>=0.1) %>% dplyr::rename(gene=feature,cluster=group, avg_log2FC=logFC)

### 测试 Marker 表达  出图
p_marker <- FeaturePlot(aging, features = c("DAZL","DDX4","MAGEA4","UTF1","FCGR3A","KIT","DMRT1","DMRTB1","STRA8","SYCP3", "SPO11", "MLH3","ZPBP","ID4","PIWIL4","UCHL1", "TNP1", "TNP2", "PRM2","SOX9", "WT1", "AMH", "PRND","FATE1","VWF","PECAM1","CDH5","DLK1", "IGF1","CYP11A1","STAR","NOTCH3","ACTA2","MYH11", "CYP26B1", "WFDC1","CD14", "CD163","C1QA","C1QC","CD8A","CD8B","PTPRC"),label = F, ncol = 5)

ggsave(
  filename = "/data/sca/user_data/4/output/pmarker.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 10000,             # 宽
  height = 9000,            # 高
  units = "px",          # 单位
  dpi = 300,              # 分辨率DPI
  plot = p_marker,
  limitsize = FALSE
)


# rda=>load rds=>read
#load("/data/sca/ref/reference.rds")
#load("/data/sca/ref/20230530-test_Lu.rda")
#load("/data/sca/ref/scRNA_zhushilabel_SE.ref2.rds")
young1 <- readRDS("/data/sca/ref/scRNA_zhushi.rds")

test_label <- t(FetchData(young1,vars = c("ident")))

test_data <- as.data.frame(GetAssayData(young1,slot = "data"))
print("test")

test_ref_list <- list(count = test_data, label = test_label)
print("test")

aging_for_SingleR <- GetAssayData(aging, slot="data") ##获取标准化矩阵
aging.hesc <- SingleR(test = aging_for_SingleR, ref = test_ref_list$count, labels = test_ref_list$label)
print("test")

aging@meta.data$labels <-aging.hesc$labels

## 将注释的label加到ident中
Idents(aging) <-aging$labels
## 定义细胞类型
aging$celltype <- aging@active.ident

plot_celltytpe <- DimPlot(aging, group.by = c("seurat_clusters", "labels"),reduction = "umap",label = T)

ggsave(
  filename = "/data/sca/user_data/4/output/annotation_result.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = plot_celltytpe,
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
  geom_bar(aes(y =freq, x= sample, fill = celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

jsonlite::toJSON(celltype_ratio, pretty = F)

ggsave(
  filename = "/data/sca/user_data/4/output/p_ratio.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = p_ratio,
  limitsize = FALSE
)


# 基因差异分析
DEG_all <- rbind.data.frame()

for(item in rev(unique(aging$celltype))){
    tmp <- subset(aging,idents = item)
    errCheck = tryCatch({
      tmp.markers <- FindMarkers(tmp,group.by = "sample", ident.1 = "expr", ident.2 = "ctrl", min.pct = 0.25)
      tmp.markers$gene <- rownames(tmp.markers)
      0
    },error = function(e){
      # print(e)
      2
    })
    if(errCheck==2){
      # print("遇到错误，结束循环")
      next
    }
    tmp.markers$cluster <- item
    DEG_all <- rbind.data.frame(DEG_all, tmp.markers)
}

write.csv(DEG_all, file = "/data/sca/user_data/4/output/different_expression_gene.csv")

p_deg_volcano <- jjVolcano(diffData = DEG_all,  tile.col = corrplot::COL2('RdBu', length(table(DEG_all$cluster))))

ggsave(
  filename = "/data/sca/user_data/4/output/p_deg_volcano.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = p_deg_volcano,
  limitsize = FALSE
)


library(clusterProfiler)

#### 针对差异基因进行通路富集分析，区分上调基因和下调基因
go.enrich=function(gene){
  eg = bitr(gene, fromType="SYMBOL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
  ego <- enrichGO(gene         = eg[,2],OrgDb= org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.1,readable = T)
  go=data.frame(ego@result)
  go$GeneRatio2<-sapply(go$GeneRatio,function(x) eval(parse(text = x)))
  return(go)
}

### 细胞类型的富集分析,针对基因的
gene.enrich=function(data){
  for (i in unique(data$cluster)){
    test=subset(data,cluster==i)
    result=go.enrich(test$gene)
    result$cluster=i
    if(i==unique(data$cluster)[1]){
      go=result
    }else{
      go=rbind(go,result)
    }
  }
  return(go)
}

### 上调基因
up=subset(DEG_all, p_val_adj<0.05 & avg_log2FC > 0.25)
up.go=gene.enrich(up)

write.csv(up.go, file = "/data/sca/user_data/4/output/upGo.csv")

upGoTopN <- up.go %>% group_by(cluster) %>% arrange(p.adjust) %>% slice_head(n = 15) %>% arrange(desc(Count))

### 下调基因
down=subset(DEG_all, p_val_adj<0.05 & avg_log2FC < -0.25)
down.go=gene.enrich(down)

write.csv(down.go, file = "/data/sca/user_data/4/output/downGo.csv")

downGoTopN <- down.go %>% group_by(cluster) %>% arrange(p.adjust) %>% slice_head(n = 15) %>% arrange(desc(Count))

# 上调基因
#纵向柱状图#
up_go_bar_plot <- ggplot(upGoTopN,aes(x=Description,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  facet_grid(ONTOLOGY~cluster, scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "Up Regulate GO Terms Enrichment")+  #设置标题
  theme_bw()

ggsave(
  filename = "/data/sca/user_data/4/output/up_go_bar_plot.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = up_go_bar_plot,
  limitsize = FALSE
)

#点图#
up_go_point_plot <- ggplot(upGoTopN,aes(x=cluster,y=reorder(Description,-pvalue),size=Count,color=-log10(pvalue)))+
  geom_point()+theme_classic()+
  theme(axis.text.x = element_text(color="black",size=13,angle=0,hjust=0.5),
        axis.text.y = element_text(color="black",size=13),
        axis.title.x = element_text( color="black",size=15),
        axis.title.y = element_text( color="black",size=15))+
  scale_color_gradient(low="lightgrey", high="red")+
  xlab("Clusters") +  #x轴标签
  ylab("Pathway") +  #y轴标签
  labs(title = "Up Regulate GO Terms Enrichment")+  #设置标题
  facet_wrap(~ONTOLOGY,ncol = 3)

ggsave(
  filename = "/data/sca/user_data/4/output/up_go_point_plot.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = up_go_point_plot,
  limitsize = FALSE
)

# 下调基因
#纵向柱状图#
down_go_bar_plot <- ggplot(downGoTopN, aes(x=Description,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  # facet_wrap("~ cluster")+
  facet_grid(ONTOLOGY~cluster, scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "Down Regulate GO Terms Enrichment")+  #设置标题
  theme_bw()

ggsave(
  filename = "/data/sca/user_data/4/output/down_go_bar_plot.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = down_go_bar_plot,
  limitsize = FALSE
)


# 下调基因
#纵向点图#
down_go_point_plot <- ggplot(downGoTopN,aes(x=cluster,y=reorder(Description,-pvalue),size=Count,color=-log10(pvalue)))+
  geom_point()+theme_classic()+
  theme(axis.text.x = element_text(color="black",size=13,angle=0,hjust=0.5),
        axis.text.y = element_text(color="black",size=13),
        axis.title.x = element_text( color="black",size=15),
        axis.title.y = element_text( color="black",size=15))+
  scale_color_gradient(low="lightgrey", high="blue")+
  xlab("Clusters") +  #x轴标签
  ylab("Pathway") +  #y轴标签
  labs(title = "Down Regulate GO Terms Enrichment")+  #设置标题
  facet_wrap(~ONTOLOGY,ncol = 3)

ggsave(
  filename = "/data/sca/user_data/4/output/down_go_point_plot.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4000,             # 宽
  height = 2000,            # 高
  units = "px",          # 单位
  dpi = 250,              # 分辨率DPI
  plot = down_go_point_plot,
  limitsize = FALSE
)


