
#-------------------------------------------------------------------------------------------------------------
#                            全部细胞绘图
#-------------------------------------------------------------------------------------------------------------
## 调用需要的R包 ##
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(stringr)
library(RColorBrewer)
library(ggsci)
library(tidydr)
library(ComplexHeatmap)
library(magrittr)
library(grid)
library(gridExtra)
library(scattermore) # 将点栅格化
library(qs)
#-------------------------------------------------------------------------------------------------------------
## 读入大群数据 ##
scRNA <- qread('/root/wangje/Project/OveryArtical/03NK_TCells/001_NKT.qs')
scRNA
# An object of class Seurat
# 21994 features across 9009 samples within 1 assay
# Active assay: RNA (21994 features, 2000 variable features)
#  2 layers present: counts, data
#  6 dimensional reductions calculated: pca, umap, Harmonypca, Harmony, umap, HarmonyUMAP3D
## 去除其中的55-60组合的数据 ##
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts,meta.data = scRNA@meta.data)
scRNA <- scRNA[,scRNA$Group1 != '55-60']
range(table(scRNA$Group1))
range(table(scRNA$sample))

###### 转换成新的anndata对象 #####
library(sceasy)
sceasy::convertFormat(scRNA,from = 'seurat',to = 'anndata',outFile = './去除55-60组合的anndata对象.h5ad')
#-------------------------------------------------------------------------------------------------------------
mark_list <- list(
  "T cells" = c('CD2','CD3G','CD3E','CD8A','CD8B','CD4'),
  "CD4+TN" = c('LEF1','TCF7','SELL','IL7R'),
  "CD4+Tem" = c("CD40LG","ANXA1","FOS","JUN"),
  "CD4+Treg" = c("FOXP3","SAT1","IL2RA","CTLA4"),
  "CD4+Tex" = c('PDCD1',"CXCL13","CD200","TNFRSF18"),
  "CD8+TN" = c('CCR7',"NELL2","CD55","KLF2"),
  "CD8+TEM" = c("GZMK","EOMES","ITM2C"),
  "CD8+Temra" = c("CX3CR1","GNLY"),
  "CD8+TEX" = c("GZMH","GZMB","LAG3","CCL4L2"),
  "NKcyto" = c("FCGR3A","FGFBP2","TYROBP"),
  "NKrest" = c('AREG',"XCL1",'KLRC1'),
  "Tys" = c('TRDV2','TRGV9','MTRNR2L8','KLRD1','TRDV1','KLRC3','CTSW'),
  "Proliferating" = c('MKI67',"STMN1","TUBA1B",'HIST1H4C'),
  TH1 = c('ZEB2',"BHLHE40",'TBX21','RUNX2','RUNX3','STAT1','TNF'),
  TH17 = c('IFNG','CCR6','CCR4','IL23R','IL17A','IL17B','FURIN','CTSH'),
  TFH = c('BCL6','EBI3','ICOS','TOX','CXCR5','STAT3','TNFSF4','SLAMF1'),
  CD8_MAIT = c('KLRB1','SLC4A10','RORC','RORA')
)

Idents(scRNA) <- scRNA$leiden_res_2.00
p2 <- DotPlot(scRNA, features = mark_list)
data <- p2$data
head(data)
data1 <- data %>% dplyr::select(avg.exp.scaled,features.plot,id)
head(data1)
data1 <- data1 %>% tidyr::pivot_wider(names_from = 'id',values_from  = 'avg.exp.scaled') %>% 
  tibble::column_to_rownames('features.plot')

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlGn")))(40) ##蓝到红
values <- seq(-2.5, 2.5, length.out = 40)[-41]
col_fun = colorRamp2(values, colors)

rowAnnotations <- rowAnnotation(
  
)

mark_list <- list(
  "T cells" = c('CD2','CD3G','CD3E','CD8A','CD8B','CD4'),
  "CD4+TN" = c('LEF1','TCF7','SELL','IL7R'),
  "CD4+Tem" = c("CD40LG","ANXA1","FOS","JUN"),
  "CD4+Treg" = c("FOXP3","SAT1","IL2RA","CTLA4"),
  "CD4+Tex" = c('PDCD1',"CXCL13","CD200","TNFRSF18"),
  "CD8+TN" = c('CCR7',"NELL2","CD55","KLF2"),
  "CD8+TEM" = c("GZMK","EOMES","ITM2C"),
  "CD8+Temra" = c("CX3CR1","GNLY"),
  "CD8+TEX" = c("GZMH","GZMB","LAG3","CCL4L2"),
  "NKcyto" = c(),
  "NKrest" = c('AREG',"XCL1",'KLRC1'),
  "Tys" = c('TRDV2','TRGV9','MTRNR2L8','KLRD1','TRDV1','KLRC3','CTSW'),
  "Proliferating" = c('MKI67',"STMN1","TUBA1B",'HIST1H4C'),
  TH1 = c('ZEB2',"BHLHE40",'TBX21','RUNX2','RUNX3','STAT1','TNF'),
  TH17 = c('IFNG','CCR6','CCR4','IL23R','IL17A','IL17B','FURIN','CTSH'),
  TFH = c('BCL6','EBI3','ICOS','TOX','CXCR5','STAT3','TNFSF4','SLAMF1'),
  CD8_MAIT = c('KLRB1','SLC4A10','RORC','RORA')
)

p2 <- DotPlot(scRNA, features = mark_list, group.by = "leiden_res_2.00")
data <- p2$data
head(data)
data1 <- data %>% dplyr::select(avg.exp.scaled,features.plot,id)
head(data1)
data1 <- data1 %>% tidyr::pivot_wider(names_from = 'id',values_from  = 'avg.exp.scaled') %>% 
  tibble::column_to_rownames('features.plot')

# 生成row标题
test <- sapply(names(mark_list),FUN = function(x){
  rep(x, length(mark_list[[x]]))
})
unlist(test)

# 绘图
data1 <- data1[as.character(0:17)]
p1 <- Heatmap(as.matrix(data1),
            col = col_fun,
            cluster_rows = F, 
            row_split = unlist(test),
            row_title_rot = 0,
            column_split = rep(LETTERS[1:ncol(data1)]),
            column_title = NULL,
            column_gap = unit(0,"cm"),
            border = TRUE,
            row_gap = unit(0,'cm'),
            cluster_columns = F)
#-------------------------------------------------------------------------------------------------------------
###### 分辨率为2.00 #####
mark_list <- list(
    "T cells" = c('CD2','CD3G','CD3E','CD8A','CD8B','CD4'),
    "NK cells" = c("FCGR3A","FGFBP2","TYROBP",'NCAM1','KLRC1')
)
p2 <- DotPlot(scRNA, features = mark_list, group.by = "leiden_res_2.00")
data <- p2$data
head(data)
data1 <- data %>% dplyr::select(avg.exp.scaled,features.plot,id)
head(data1)
data1 <- data1 %>% tidyr::pivot_wider(names_from = 'id',values_from  = 'avg.exp.scaled') %>% 
  tibble::column_to_rownames('features.plot')

###### 分辨率为1.00 #####
mark_list <- list(
    "T cells" = c('CD2','CD3G','CD3E','CD8A','CD8B','CD4'),
    "NK cells" = c("FCGR3A","FGFBP2","TYROBP",'NCAM1','KLRC1')
)
p2 <- DotPlot(scRNA, features = mark_list, group.by = "leiden_res_1.00")
data <- p2$data
head(data)
data1 <- data %>% dplyr::select(avg.exp.scaled,features.plot,id)
head(data1)
data1 <- data1 %>% tidyr::pivot_wider(names_from = 'id',values_from  = 'avg.exp.scaled') %>% 
  tibble::column_to_rownames('features.plot')
data1 <- data1[as.character(0:max(as.numeric(colnames(data1))))]
pheatmap(data1,cluster_cols = FALSE)

#-------------------------------------------------------------------------------------------------------------
###### 细胞命名 #####
scRNA$celltype <- ifelse(scRNA$leiden_res_1.00 %in% c(2,0,7,8),'CD4+ T cells',ifelse(scRNA$leiden_res_1.00 %in% c(6),'CD56+CD16-NK','CD8+ T cells'))
scRNA$celltype <- ifelse(scRNA$leiden_res_2.00 %in% c(7),'CD56+CD16-NK',
    ifelse(scRNA$leiden_res_2.00 %in% c(8),'CD56+CD16+NK',
    ifelse(scRNA$leiden_res_2.00 %in% c(13),'Treg',scRNA$celltype)))

#-------------------------------------------------------------------------------------------------------------
## 添加分组信息 ##
scRNA$sample <- stringr::str_split_fixed(Cells(scRNA), "_[A|T|G|C].*", n = 2)[, 1]
scRNA@meta.data %<>% dplyr::mutate(
    Group1 = dplyr::case_when(
        sample %in% c("young1", "young3", "young4", "na1111", "na1122", "na-426", "na-1010") ~ "18-35",
        # sample %in% c('na1122','HRS421453','na-426','na-1010') ~ '32-35',
        sample %in% c("middle2", "middle3", "middle4") ~ "37-39",
        sample %in% c("na1128", "na1117", "HRS421449", "HRS421450") ~ "39-44",
        sample %in% c("old1", "old2", "old3") ~ "47-49",
        sample %in% c("na1129", "NA_728", "na-412") ~ "55-60"
    ),
    Group1 = factor(Group1, levels = c("18-35", "37-39", "39-44", "47-49", "55-60")),
    Group2 = dplyr::case_when(
        sample %in% c("young1", "young3", "young4", "na1111") ~ "18-29",
        sample %in% c("na1122", "na-426", "na-1010") ~ "32-35",
        sample %in% c("middle2", "middle3", "middle4") ~ "37-39",
        sample %in% c("na1128", "HRS421449", "HRS421450", "na1117") ~ "39-44",
        sample %in% c("old1", "old2", "old3") ~ "47-49",
        sample %in% c("na1129", "NA_728", "na-412") ~ "55-60"
    ),
    Group2 = factor(Group2, levels = c("18-29", "32-35", "37-39", "39-44", "47-49", "55-60")),
    Menopause = dplyr::case_when(
        sample %in% c("na1129", "NA_728", "na-412") ~ "Yes",
        TRUE ~ "No"
    ),
    Menopause = factor(Menopause, levels = c("No", "Yes"))
)
# 添加Group4分组
if(is.factor(scRNA$Group1)) scRNA$Group1 <- varhandle::unfactor(scRNA$Group1)
scRNA$Group4 <- ifelse(scRNA$Group1 %in% c('37-39','39-44'),'37-44',scRNA$Group1)

cat('分组结果')
print(unique(scRNA$sample))
print(unique(scRNA$Group1))
print(unique(scRNA$Group4))    

#-------------------------------------------------------------------------------------------------------------
if(!is.factor(scRNA$celltype)) scRNA$celltype <- factor(scRNA$celltype, levels = names(sort(table(scRNA$celltype),decreasing = TRUE)))
print(levels(scRNA$celltype))
## 输出样本细胞数 ##
df <- scRNA@meta.data %>% group_by(sample,Group1,Group4) %>% dplyr::summarise(cells = n()) %>%dplyr::arrange(Group1,Group4)
p <- ggtexttable(df, rows = NULL)
ggsave(p, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/样本细胞数.png", width = 6, height = 10,bg = 'white')

## 查看数据的分布情况 ##
scRNA$new_group <- paste0(scRNA$sample,'_',scRNA$Group1)
if(!("pecent.mt" %in% colnames(scRNA@meta.data))| is.null(scRNA$percent.mt)) scRNA$percent <- PercentageFeatureSet(scRNA, pattern = "^MT-")
p <- VlnPlot(scRNA, group.by = 'new_group', features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3,pt.size = 0)
ggsave(p, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/样本细胞数分布图.png", width = 18, height = 5,bg = 'white')

#-------------------------------------------------------------------------------------------------------------
DimPlot2 <- function(srt,
                     group.by,
                     pt.size = 0.01,
                     reduction = "umap",
                     xlab,
                     ylab,
                     startNum = 0,
                     # 启始序号
                     label.size = 4,
                     title = NULL,
                     raster = FALSE,
                     split.by = NULL,
                     ncol = NULL,
                     # 分面展示列数
                     cols = c(
                       RColorBrewer::brewer.pal(n = 7, name = "Dark2"),
                       RColorBrewer::brewer.pal(n = 7, name = "Accent")
                     ),
                     pelate = NULL,
                     use_theme = "theme_dr",
                     ...) {
  draw_number_circle <- function(data, params, size) {
    grobTree(
      pointsGrob(
        x = 0.5,
        y = 0.5,
        size = unit(1.6, "char"),
        pch = 16,
        gp = gpar(
          col = alpha(data$colour %||% "grey50", data$alpha),
          fill = alpha(data$fill %||% "grey50", data$alpha),
          lwd = (data$linewidth %||% 0.5) * .pt,
          lty = data$linetype %||% 1
        )
      ),
      textGrob(
        label = data$label,
        x = rep(0.5, 3),
        y = rep(0.5, 3),
        gp = gpar(col = "black")
      )
    )
  }
  # 查看分组类型是否在meta数据中,如果不是因子水平怎添加因子水平
  if (!(group.by %in% names(srt@meta.data)))
    stop(sprintf("%s is not in Seurat meta.data", group.by))
  if (!(is.factor(srt@meta.data[[group.by]])))
    srt@meta.data[[group.by]] <- factor(srt@meta.data[[group.by]], levels = names(sort(table(srt@meta.data[[group.by]]), decreasing = TRUE)))
  # 获取绘图的坐标数据
  embedding <-
    Seurat::Embeddings(srt, reduction = reduction) %>% as.data.frame()
  colnames(embedding) <- c('UMAP_1', 'UMAP_2')
  meta.data <- srt@meta.data
  data <- cbind(embedding, meta.data)
  # 添加序号
  # id 根据level水平进行添加
  # print(levels(data[[group.by]]))
  if (startNum == 0) {
    tt = data.frame(levels(data[[group.by]]), startNum:(length(unique(data[[group.by]])) - 1))
  } else if (startNum == 1) {
    tt = data.frame(levels(data[[group.by]]), startNum:(length(unique(data[[group.by]]))))
  } else{
    warning("The start number can only be 0 or 1")
    startNum <- 0
    tt = data.frame(levels(data[[group.by]]), startNum:(length(unique(data[[group.by]])) - 1))
  }
  colnames(tt) <- c(group.by, 'id')
  tt[[group.by]] <-
    factor(tt[[group.by]], levels = levels(data[[group.by]]))
  data <- dplyr::left_join(data, tt, by = group.by)
  # print(head(data))
  # 计算label添加的位置
  label_id <- data %>% group_by(id) %>%
    dplyr::summarise(x_m = median(UMAP_1), y_m = median(UMAP_2)) # 添加在每个分区中的中心位置
  label_id$id <- factor(label_id$id, levels = unique(data$id))
  # 进行绘图
  p <- NULL
  if (is.null(xlab))
    xlab <- "UMAP1"
  if (is.null(ylab))
    ylab <- "UMAP"
  if (!raster) {
    p <- ggplot(data = data, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = !!sym(group.by)), key_glyph = draw_number_circle, size = pt.size) +
      geom_label(
        data = label_id,
        aes(x = x_m, y = y_m, label = id),
        size = label.size,
        alpha = 0,
        label.r = unit(0.25, "cm"),
        label.size = NA
      ) +
      guides(color = guide_legend(override.aes = list(label = label_id$id))) +
      labs(x = xlab,
           y = ylab,
           title = paste0("n=", format(
             dim(srt)[2], big.mark = ",", scientific = FALSE
           )))
  } else{
    p <- ggplot(data = data, aes(x = UMAP_1, y = UMAP_2)) +
      geom_scattermore(aes(color = !!sym(group.by)), key_glyph = draw_number_circle, size = pt.size) +
      geom_label(
        data = label_id,
        aes(x = x_m, y = y_m, label = id),
        size = label.size,
        alpha = 0,
        label.r = unit(0.25, "cm"),
        label.size = NA
      ) +
      guides(color = guide_legend(override.aes = list(label = label_id$id))) +
      labs(x = xlab,
           y = ylab,
           title = paste0("n=", format(
             dim(srt)[2], big.mark = ",", scientific = FALSE
           )))
  }
  # 是否进行分面
  if (!(is.null(split.by))) {
    if (!(split.by %in% names(data)))
      stop(sprintf("%s is not in data, please check!", split.by))
    p <-
      p + facet_wrap(facets = vars(!!sym(x = split.by)),
                     ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                       length(x = unique(x = data[, split.by]))
                     },
                     scales = "free") +
      labs(title = NULL)
  }
  
  # 添加theme
  p <- p & theme_classic() &
    theme(
      legend.title = element_blank(),
      panel.grid = element_blank(),
      legend.text = element_text(size = label.size + 7, color = "black"),
      axis.title = element_text(size = label.size + 8, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = label.size + 6, colour = "black"),
      axis.text = element_text(size = label.size + 6, colour = "black"),
      ...
    )
  if (!(is.null(cols)))
    p <- p + ggplot2::scale_color_manual(values = cols)
  if (is.null(cols) &&
      !(is.null(palette)))
    p <- p + ggplot2::scale_color_brewer(palette = palette)
  return(p)
}

###### 绘制DimPlot #####
scRNA$Group1 <- factor(scRNA$Group1, levels = c("18-35", "37-39", "39-44", "47-49"))
scRNA$Group4 <- factor(scRNA$Group4, levels = c("18-35", "37-44", "47-49"))
## celltype ##
cols <- c('#9e1a42','#d44351','#ff79a1',"#7050b1",'#b185df','#85d1f8','#386da9','#133d56','#245801','#528d1f','#a4d874')
if(!("celltype" %in% colnames(scRNA@meta.data))) stop('no celltype column in Seurat meta.data')
p1 <- DimPlot2(scRNA,group.by = "celltype",reduction = "umap",startNum = 0,xlab = "UMAP1",ylab = "UMAP2",split.by = NULL,ncol = NULL,label.size = 4,cols = RColorBrewer::brewer.pal(9,'Set1'))
p1 <- p1 + theme_classic()+
  theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),legend.title = element_blank())
# 保存图片
if(!dir.exists('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群'))dir.create('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群',recursive = TRUE)
ggsave(p1, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot2_celltype.png", width = 6, height = 4,bg = 'white')
ggsave(p1, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot2_celltype.pdf", width = 6, height = 4,bg = 'white')

## sample ##
p2 <- DimPlot(scRNA, group.by = 'sample',cols = c(RColorBrewer::brewer.pal(n = 8,name = 'Set1'),RColorBrewer::brewer.pal(n = 12,name = 'Paired')), reduction = "umap")+
    labs(x='UMAP1',y='UMAP2')+ guides(color = guide_legend(ncol=2,override.aes = list(size = 5))) + theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
    theme(legend.title = element_blank(),panel.grid = element_blank())
ggsave(p2, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_sample.png", width = 6, height = 4,bg = 'white')
ggsave(p2, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_sample.pdf", width = 6, height = 4,bg = 'white')
        
## Group ##   
print(unique(scRNA$Group1))   #[1] "39-44" "18-35" "37-39" "47-49"
p3 <- DimPlot(scRNA, group.by = 'Group1',cols = c('#dfb67d','#cb8eab','#869746','#4f8799'), reduction = "umap")+
        labs(x='UMAP1',y='UMAP2')+ guides(color = guide_legend(ncol=1,override.aes = list(size = 5))) + theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
    theme(legend.title = element_blank(),panel.grid = element_blank())
ggsave(p3, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_Group1.png", width = 5, height = 4,bg = 'white')
ggsave(p3, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_Group1.pdf", width = 5, height = 4,bg = 'white')

#-------------------------------------------------------------------------------------------------------------
###### 绘制split dimplot #####
library(SCP)
# 绘制不同分组的celltype
p4 <- SCP::CellDimPlot(scRNA, group.by = 'celltype',reduction = 'umap',split.by = 'Group1',palette = 'Set1',xlab = 'UMAP1',ylab = 'UMAP2',theme_use = 'theme_blank',bg_color = 'white',ncol = 4)
ggsave(p4, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_Group1_split.png", width = 24, height = 4,bg = 'white')
ggsave(p4, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_Group1_split.pdf", width = 24, height = 4,bg = 'white')

p4 <- SCP::CellDimPlot(scRNA, group.by = 'celltype',reduction = 'umap',split.by = 'Group4',palette = 'Set1',xlab = 'UMAP1',ylab = 'UMAP2',theme_use = 'theme_blank',bg_color = 'white',ncol = 3)
ggsave(p4, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_Group4_split.png", width = 18, height = 4,bg = 'white')
ggsave(p4, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_Group4_split.pdf", width = 18, height = 4,bg = 'white')

# 绘制每个样本的celltype split
print(unique(scRNA$sample))
print(length(unique(scRNA$sample)))
p5 <- SCP::CellDimPlot(scRNA, group.by = 'celltype',reduction = 'umap',split.by = 'sample',palette = 'Set1',xlab = 'UMAP1',ylab = 'UMAP2',theme_use = 'theme_blank',bg_color = 'white',ncol = 4)
ggsave(p5, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_sample_split.png", width = 30, height = 12,bg = 'white', limitsize = FALSE)
ggsave(p5, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/DimPlot_sample_split.pdf", width = 30, height = 12,bg = 'white',limitsize = FALSE)
#-------------------------------------------------------------------------------------------------------------
###### 绘制热图 #####
## 计算差异基因 ##
if(!is.factor(scRNA$celltype)) scRNA$celltype <- factor(scRNA$celltype,levels = names(sort(table(scRNA$celltype))))
Idents(scRNA) <- scRNA$celltype
scRNA.FindAllMarks <- FindAllMarkers(scRNA,only.pos = TRUE, min.pct = 0.25,pval = 0.05,logfc.threshold = 0.25)
# 保留获得的差异基因结果
write.table(scRNA.FindAllMarks,file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/FindAllMarkers_onlyPos.txt",quote = FALSE,sep = '\t')
# 过滤差异基因
# 保存Top10基因的结果
Top10.genes <- scRNA.FindAllMarks %>% group_by(cluster) %>% arrange(desc(avg_log2FC), p_val_adj) %>% dplyr::slice(1:10)
write.table(Top10.genes,file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/Top10_genes.txt",quote = FALSE,sep = '\t')
## 绘制热图 ##
scRNA <- ScaleData(scRNA)
p6 <- DoHeatmap(scRNA,features = Top10.genes$gene)+scale_fill_viridis_c()
ggsave(p6, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/NKT大群TOP10基因热图.png", width = 7, height = 8,bg = 'white', limitsize = FALSE)

#-------------------------------------------------------------------------------------------------------------
###### 绘制marker基因的FeaturePlot #####
plotMarkerFeatures <- function(srt, mark_list,reduction,
  cols = c("#FFEFD5", "#E6E6FA", "#87CEFA", "#6495ED", "#4169E1", "#0000CD", "#000080")) {  
  plist <- list()  
  if (!is.list(mark_list)) {  
    stop("mark_list must be a list.")  
  }  
  
  # 确保srt对象包含必要的属性  
  if (!(inherits(srt,'Seurat'))) {  
    stop("srt must be a data frame with 'rownames' and 'umap' attributes.")  
  }  
  
  for (i in seq_along(mark_list)) {  
    if (length(mark_list[[i]]) == 0) {  
      warning(paste("Sublist", i, "of mark_list is empty."))  
      next  
    }  
    
    for (j in mark_list[[i]]) {  
      if (j %in% rownames(srt)) {  
        p <- FeaturePlot(srt,  # 假设FeaturePlot是某个包中的函数  
                         features = j,  
                         reduction = reduction,  
                         max.cutoff = 1.5,  
                         raster = FALSE,  
                         pt.size = 0.03,  
                         cols = cols  
        ) +  
          scale_x_continuous("") +  
          scale_y_continuous("") +  
          theme_bw() +  
          theme(  
            panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank(),  
            axis.ticks = element_blank(),  
            axis.text = element_blank(),  
            legend.position = "none",  
            plot.title = element_text(hjust = 0.5, size = 18)  
          ) +  
          ggtitle(paste0(names(mark_list)[i], " (", j, ")"))  # 修改此处  
        plist[[paste0(i, "|", j)]] <- p  
      }  
    }  
  }  
  return(plist)  
}

mark_list <- list(
    "T cells" = c('CD3G','CD3E','CD4','CD8A','CD8B'),
    'NK cells' = c("FCGR3A","FGFBP2","TYROBP",'NCAM1','KLRC1')
)
# cols <- (RColorBrewer::brewer.pal(n=8,name = 'Blues'))
plist <- plotMarkerFeatures(scRNA, mark_list = mark_list, reduction = 'umap')
ggsave(patchwork::wrap_plots(plist,ncol = 5), filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞marker基因FeaturePlot图.png", width = 15, height = 6,bg = 'white', limitsize = FALSE)
# ggsave(patchwork::wrap_plots(plist,ncol = 6), filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞marker基因FeaturePlot图.pdf", width = 10, height = 4,bg = 'white',limitsize = FALSE)

#-------------------------------------------------------------------------------------------------------------
###### 绘制细胞比例变化箱线图 #####
library(ggpubr)
library(patchwork)
library(Seurat)
library(dplyr)
library(rstatix)
library(purrr)
library(varhandle)

## 计算细胞变化比例 ##
prop_caculate <- function(srt,celltype,group1,group2){
  srt@meta.data %>% dplyr::group_by(.data[[celltype]],.data[[group1]],.data[[group2]]) %>% 
    dplyr::summarise(n = n()) %>% dplyr::left_join(
      srt@meta.data  %>% dplyr::group_by(.data[[group1]]) %>% dplyr::summarise(total = n()),by=group1
    )%>% dplyr::mutate(
      prop = n/total
    )
}
if(interactive()){
  df <- prop_caculate(scRNA,celltype = 'celltype',group1 = "sample",group2 = "Group1")}
# 写出文件
write.csv(df,'/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞比例变化.csv',row.names = FALSE, quote=FALSE)
# 去除其中的卵细胞和颗粒细胞
df <- df[!(df$celltype %in% c('Oocyte','Granulosa')),]
if(is.factor(df$celltype)) df$celltype <- varhandle::unfactor(df$celltype)

## 绘制比例变化箱线图 ##
plot_boxplot <- function(df, x, y, fill, compair, hide.ns = TRUE) {
  result <- map(unique(df[[x]]), function(i) {
    print(i)
    df_sub <- df %>% filter(!!sym(x) == i)
    # 计算比例
    stat.test <- df_sub %>%
      group_by(!!sym(x)) %>%
      wilcox_test(as.formula(paste0(y,'~',fill)), comparisons = compair) %>%
      add_xy_position(x = x, dodge = 0.8)
    
    if (hide.ns) {
      stat.test$new_signif <- case_when(
        0 <= stat.test$p & stat.test$p < 0.01 ~ "***",
        0.01 <= stat.test$p & stat.test$p < 0.05 ~ "**",
        0.05 <= stat.test$p & stat.test$p < 0.1 ~ "*",
        TRUE ~ "ns"
      )
      stat.test <- stat.test %>% filter(new_signif != "ns")
    }
    
    p2 <- ggplot(df_sub) + 
      geom_boxplot(aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]]),
                   position = position_dodge(width = 0.8),
                   outlier.shape = NA, color = "black") + 
      geom_jitter(aes(x = !!sym(x), y = !!sym(y), fill = !!sym(fill)),
                  color = "black", position = position_dodge(width = 0.8), pch = 21, size = 3)+
      ggpubr::stat_kruskal_test(aes(x = .data[[x]], y = .data[[y]],group=.data[[fill]]), label.x = 0.8,label.y=max(df_sub$prop)+0.1,
                                label = "p={p.format}")+
      stat_pvalue_manual(stat.test, label = "new_signif", tip.length = 0.00, size = 8, hide.ns = FALSE) +
      theme_classic(base_size = 20, base_line_size = 0.6) +
      # scale_fill_manual(values = c(brewer.pal(n = 12, name = "Paired"))) +
      scale_fill_manual(values = c('#dfb67d','#cb8eab','#869746','#4f8799'))+
      labs(x="",y="Fraction") +
      ggtitle(i)+
      theme(
        legend.key.size = unit(1, 'cm'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 17, colour = "black"),
        # axis.text.x = element_text(size = 17, colour = "black",hjust = 1,vjust = 1,angle = 35),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 15, colour = "black"),
        axis.ticks.y = element_line(linewidth = 0.8),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5, colour = "black", size = 18, angle = 0)
      ) +
      scale_x_discrete(expand = c(0.42, 0))
    return(p2)
  })
  return(result)
}

# 运行函数
compair = combn(unique(if(is.factor(scRNA$Group1)) unique(varhandle::unfactor(scRNA$Group1)) else unique(scRNA$Group1)),m=2,simplify = FALSE) # 生长两两比较组合的列表
plist <- plot_boxplot(df, x="celltype",y="prop",fill="Group1",compair = compair)
# 保存图片
ggsave(ggpubr::ggarrange(plotlist = plist,ncol = 5,nrow = 1,common.legend = TRUE, legend = 'bottom'), filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞比例箱线图.png", width = 15, height = 4,bg = 'white', limitsize = FALSE)
ggsave(ggpubr::ggarrange(plotlist = plist,ncol = 5,nrow = 1,common.legend = TRUE, legend = 'bottom'), filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞比例箱线图.pdf", width = 15, height = 4,bg = 'white', limitsize = FALSE)

#-------------------------------------------------------------------------------------------------------------
###### 绘制3个分组的细胞比例箱线图 #####
if(interactive()){
  df <- prop_caculate(scRNA,celltype = 'celltype',group1 = "sample",group2 = "Group4")}
# 写出文件
write.csv(df,'/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞比例变化Group2.csv',row.names = FALSE, quote=FALSE)
# 去除其中的卵细胞和颗粒细胞
df <- df[!(df$celltype %in% c('Oocyte','Granulosa')),]
if(is.factor(df$celltype)) df$celltype <- varhandle::unfactor(df$celltype)
# 运行函数
compair = combn(unique(if(is.factor(scRNA$Group4)) unique(varhandle::unfactor(scRNA$Group4)) else unique(scRNA$Group4)),m=2,simplify = FALSE) # 生长两两比较组合的列表
plist <- plot_boxplot(df, x="celltype",y="prop",fill="Group4",compair = compair)
# 保存图片
ggsave(ggpubr::ggarrange(plotlist = plist,ncol = 5,nrow = 1,common.legend = TRUE, legend = 'bottom'), filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞比例箱线图Group2.png", width = 13, height = 4,bg = 'white', limitsize = FALSE)
ggsave(ggpubr::ggarrange(plotlist = plist,ncol = 5,nrow = 1,common.legend = TRUE, legend = 'bottom'), filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞比例箱线图Group2.pdf", width = 13, height = 4,bg = 'white', limitsize = FALSE)
#-------------------------------------------------------------------------------------------------------------
###### 计算并绘制Ro/e或者OR值热图 #####
## 计算OR值 ##
library(Seurat)
library(dplyr)
library(rlang)
library(plyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(data.table)

#  编写函数
#**************************************************************************
Caculate.OR <- function(srt, 
                        loc, 
                        celltype){
    if(!inherits(srt, "Seurat")){
        stop("input file is not Seurat Object, Please check!")
    }
    # 计算不同组织中不同细胞类型的细胞量
    count.dist.melt.tb <-             as.data.frame(table(srt@meta.data[[celltype]],srt@meta.data[[loc]]))
    colnames(count.dist.melt.tb) <- c('rid','cid','count')
    count.dist.melt.tb$rid <- unfactor(count.dist.melt.tb$rid)
    count.dist.melt.tb$cid <- unfactor(count.dist.melt.tb$cid)
    sum.col <- table(srt@meta.data[[loc]])
    sum.row <- table(srt@meta.data[[celltype]])
    count.dist.melt.ext.tb <-
     as.data.table(ldply(seq_len(nrow(
            count.dist.melt.tb
        )), function(i) {
            this.row <- count.dist.melt.tb$rid[i]
            this.col <- count.dist.melt.tb$cid[i]
            this.c <- count.dist.melt.tb$count[i]
            other.col.c <- sum.col[this.col] - this.c
            this.m <- matrix(c(
                this.c,
                sum.row[this.row] - this.c,
                other.col.c,
                sum(sum.col) - sum.row[this.row] - other.col.c
            ),
            ncol = 2)
            res.test <- fisher.test(this.m)
            data.frame(
                rid = this.row,
                cid = this.col,
                p.value = res.test$p.value,
                OR = res.test$estimate
            )
        }))
    count.dist.melt.ext.tb[, adj.p.value := p.adjust(p.value, "BH")]
    return(as.data.frame(count.dist.melt.ext.tb))
}

# 计算OR值
or <- Caculate.OR(scRNA,loc = "Group1",celltype = "celltype")
or2 <- Caculate.OR(scRNA,loc = "Group4",celltype = "celltype")
# 写出文件
write.table(or,file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group1_OR.txt",quote = F,sep = "\t",row.names = F)
write.table(or2,file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group4_OR.txt",quote = F,sep = "\t",row.names = F)
###### 绘制OR值热图 #####
Plot_OR <- function(inputfile){
    test <- pivot_wider(inputfile,
        names_from = 'cid',  # 指定哪些列的值应该变成新列的名称
        values_from = 'OR', # 指定哪些列的值应该填充到新列中
        id_cols = 'rid'         # 指定哪些列应该保持不变，作为新数据框的行标识
      ) -> test2
    test <- as.data.frame(test)
    test2 = test2 %>% tibble::column_to_rownames(var = 'rid')
    p1 <- pheatmap(test2, 
            scale = "row", # 按行归一化，查看因子在不同样本中的分布情况
            cluster_cols = FALSE, clustering_distance_rows = "correlation", #取消列聚类，表示行聚类使用皮尔森相关系数聚类
            treeheight_row = 30, # 设置行聚类树高
            cutree_rows =3, #根据样品列聚类情况将热图的行方向隔开为3份
            cellwidth = 30,cellheight = 15, # 设置热图方块宽度和高度
            fontsize_number = 8, #热图上数值的字体大小
            number_color="blue", #热图上数值的字体颜色
            display_numbers = TRUE,
            main="OR(odds ratios)", # 设置图形标题
            show_colnames = T, # 设置行列标签的显示
            show_rownames = T,
            border="white", 
            legend = T, # FALSE去除图例; T显示图例
            # legend_breaks=c(0,0.5,1,1.5,2,2.5,3), # 设置图例的范围
            fontsize_row = 10, # 分别设置行列标签字体大小
            fontsize_col = 10,
            angle_col = "45", # 设置标签显示角度
         )
    return(p1)
}
# 绘图并保存文件
png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group1_OR值热图.png', height = 6,width = 5, unit='in', res=300) 
Plot_OR(or)
dev.off()

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group1_OR值热图.pdf', height = 6,width = 5, unit='in', res=300) 
Plot_OR(or)
dev.off()

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group4_OR值热图.png', height = 6,width = 5, unit='in', res=300) 
Plot_OR(or2)
dev.off()

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group4_OR值热图.pdf', height = 6,width = 5, unit='in', res=300) 
Plot_OR(or2)
dev.off()

#-------------------------------------------------------------------------------------------------------------
###### 绘制Ro/e图 #####
# 计算Ro/e值
Caculate.Roe <- function(srt,loc,celltype){
    mtx <- table(srt@meta.data[[celltype]], srt@meta.data[[loc]])
    # 计算卡方值
    chisq <- chisq.test(mtx)
    Roe <- chisq$observed/chisq$expected
    return(Roe)
}

# 计算Ro/e值
Roe <- Caculate.Roe(scRNA,  loc = "Group1", celltype = "celltype")
Roe2 <- Caculate.Roe(scRNA,  loc = "Group4", celltype = "celltype")
write.table(as.data.frame(Roe), file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group1_Roe值.txt" , sep = '\t', quote = F, row.names = F)
write.table(as.data.frame(Roe2), file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group4_Roe值.txt" , sep = '\t', quote = F, row.names = F)

# 绘制Ro/e热图
# Roe > 1, "+++"
# 1 <= Roe & Roe >0.8 ,'++'
# 0.8 <= Roe &  Roe >= 0.2,'+'
# Roe > 0 & Roe < 0.2, '+/-'
# 其他：'-'
# pheatmap
p4 <- pheatmap(Roe, 
            scale = "row", # 按行归一化，查看因子在不同样本中的分布情况
            cluster_cols = FALSE, clustering_distance_rows = "correlation", #取消列聚类，表示行聚类使用皮尔森相关系数聚类
            treeheight_row = 30, # 设置行聚类树高
            cutree_rows =3, #根据样品列聚类情况将热图的行方向隔开为3份
            cellwidth = 30,cellheight = 15, # 设置热图方块宽度和高度
            fontsize_number = 8, #热图上数值的字体大小
            number_color="blue", #热图上数值的字体颜色
            # number_format="%.1", #热图上数值的字体类型
            display_numbers = matrix(ifelse(Roe > 1, "+++", ifelse(1 <= Roe & Roe >0.8 ,'++',
                ifelse(0.8 <= Roe &  Roe >= 0.2,'+',
                ifelse(Roe > 0 & Roe < 0.2, '+/-','-'
            )))), nrow(Roe)), #设置热图区分标记
            main="Ro/e", # 设置图形标题
            show_colnames = T, # 设置行列标签的显示
            show_rownames = T,
            border="white", # 设置边框为白色
            legend = T, # FALSE去除图例; T显示图例
            # legend_breaks=c(0,0.5,1,1.5,2,2.5,3), # 设置图例的范围
            fontsize_row = 10, # 分别设置行列标签字体大小
            fontsize_col = 10,
            angle_col = "45", # 设置标签显示角度
         )
p4 <- ggplotify::as.ggplot(p4)
png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group1_Roe值热图.png', height = 6,width = 7, unit='in', res=300) 
print(p4)
dev.off()

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group1_Roe值热图.pdf', height = 6,width = 7, unit='in', res=300) 
print(p4)
dev.off()

p4 <- pheatmap(as.matrix(Roe2), 
            scale = "row", # 按行归一化，查看因子在不同样本中的分布情况
            cluster_cols = FALSE, clustering_distance_rows = "correlation", #取消列聚类，表示行聚类使用皮尔森相关系数聚类
            treeheight_row = 30, # 设置行聚类树高
            cutree_rows =3, #根据样品列聚类情况将热图的行方向隔开为3份
            cellwidth = 30,cellheight = 15, # 设置热图方块宽度和高度
            fontsize_number = 8, #热图上数值的字体大小
            number_color="blue", #热图上数值的字体颜色
            # number_format="%.1", #热图上数值的字体类型
            display_numbers = matrix(ifelse(Roe2 > 1, "+++", ifelse(1 <= Roe2 & Roe2 >0.8 ,'++',
                ifelse(0.8 <= Roe2 &  Roe2 >= 0.2,'+',
                ifelse(Roe2 > 0 & Roe2 < 0.2, '+/-','-'
            )))), nrow(Roe2)), #设置热图区分标记
            main="Ro/e", # 设置图形标题
            show_colnames = T, # 设置行列标签的显示
            show_rownames = T,
            border="white", # 设置边框为白色
            legend = T, # FALSE去除图例; T显示图例
            # legend_breaks=c(0,0.5,1,1.5,2,2.5,3), # 设置图例的范围
            fontsize_row = 10, # 分别设置行列标签字体大小
            fontsize_col = 10,
            angle_col = "45", # 设置标签显示角度
         )
p4 <- ggplotify::as.ggplot(p4)

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group4_Roe值热图.png', height = 6,width = 7, unit='in', res=300) 
print(p4)
dev.off()

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype_Group4_Roe值热图.pdf', height = 6,width = 7, unit='in', res=300) 
print(p4)
dev.off()

#-------------------------------------------------------------------------------------------------------------
###### 绘制百分比柱形图 #####
## celltype细胞数 ##
df <- scRNA@meta.data %>% group_by(celltype) %>% dplyr::summarise(cells = n())
p <- ggtexttable(df, rows = NULL)
ggsave(p, filename = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/celltype样本细胞数.png", width = 4, height = 4,bg = 'white')


library(dplyr)
library(ggplot2)
library(Seurat)
library(ggalluvial)
library(RColorBrewer)
library(ggh4x)
library(qs)

## 计算比例 ##
barplot_fraction <- function(srt, celltype_col, group_col, sample_col) {  
  if (inherits(srt, "Seurat")) {  
    counts <- srt@meta.data %>%  
      dplyr::group_by(.data[[group_col]], .data[[sample_col]], .data[[celltype_col]]) %>%  
      dplyr::summarise(count = n(), .groups = "drop") %>% 
      dplyr::left_join(srt@meta.data %>% dplyr::group_by(.data[[sample_col]]) %>% dplyr::count(name = "total"),
                       by = sample_col) %>% dplyr::mutate(Freq = count/total)
    return(counts)
  } else {  
    stop("srt 不是一个 Seurat 对象")  
  }  
}  

## 绘图 ##
plot_alluvial_barplot <- function(data,
                                  x, 
                                  y, 
                                  stratum_var, 
                                  alluvium_var,  
                                  fill_var, 
                                  fill_colors = NULL, 
                                  strip_fill,
                                  ncolor = 8, 
                                  palette = "Dark2", 
                                  split.by = NULL, 
                                  label.size=16,
                                  use_ggh4x = FALSE, ...) {  
  # 检查fill_colors参数  
  if (is.null(fill_colors)) {  
    unique_fills <- unique(data[[fill_var]])  
    n_unique <- length(unique_fills)  
    fill_colors <- RColorBrewer::brewer.pal(n = min(n_unique, 12), name = "Set1")[1:n_unique]  
    if (n_unique > length(fill_colors)) {  
      warning("Not enough colors provided, reusing colors.")  
    }  
  }  
  # 初始化绘图  
  pp <- ggplot(data = data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill_var]],  
                                stratum = .data[[stratum_var]], alluvium = .data[[alluvium_var]])) +  
    geom_flow(stat = "alluvium", lode.guidance = "frontback", curve_type = "linear", 
              alpha = 0.5, width = 0.7, color = "white") +  
    geom_stratum(alpha = 1, color = "black", width = 0.7)  
  # 是否进行split  
  if (!is.null(split.by) && use_ggh4x) {  
    require(ggh4x)  # 确保已加载ggh4x包  
    class_strip <- strip_themed(background_x = elem_list_rect(fill = strip_fill)) 
    pp <- pp + ggh4x::facet_grid2(cols = vars(!!sym(split.by)), scales = "free", space = "free_x", switch = "y", strip = class_strip)  
  }
  # 修改主题
  pp <- pp + 
    labs(x = "", y = "% Fraction") +
    theme_bw()+
    theme(
      panel.background = element_blank(),
      strip.background = element_rect(colour = "black"),
      panel.grid = element_blank(),
      legend.key.size = unit(0.6,'cm'),
      panel.border = element_rect(colour = "black", linewidth = 1.3),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 1, vjust = 0.5, size = label.size+2, color = "black"),
      # panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid"),
      legend.text = element_text(size = label.size, color = "black"),
      axis.title.y = element_text(size = label.size, colour = "black"),
      axis.text.y = element_text(size = label.size, colour = "black"),
      strip.text = element_text(size=label.size+2,color="black"),
      legend.position = 'right',
      axis.text.x = element_text(size = label.size, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
      # axis.text.x = element_text(size = 20, colour = "black", angle = 45, hjust = 1, vjust = 1),
      axis.ticks.y = element_line(linewidth = 1.3),
      axis.ticks.x = element_line(linewidth = 1.3),
      axis.ticks.length.y = unit(0.4, "cm"),
      axis.ticks.length.x = unit(0.4, "cm"),
      # axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      ...
    ) +
    scale_fill_manual(values=fill_colors)
  return(pp)  
}
## 运行函数
### 样本与celltype
### Group1 4个分组
prop <- barplot_fraction(srt = scRNA,
                 celltype_col = "celltype",
                 group_col = "Group1",
                 sample_col = "sample")
if(is.numeric(prop$Freq)) prop$Freq <- 100 * prop$Freq else prop$Freq <- as.numeric(prop$Freq) * 100               
head(prop)
# 保存文件
write.table(prop, file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图sample_celltype.txt", sep = "\t", quote = F, row.names = F)

p1 <- plot_alluvial_barplot(data = prop, 
                      x = 'sample', 
                      y = 'Freq', 
                      fill_var = 'celltype',
                      fill_colors = RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                      strip_fill = c('#dfb67d','#cb8eab','#869746','#4f8799'),                                      
                      stratum_var = 'celltype', 
                      alluvium_var  ='celltype',
                      split.by = "Group1",
                      use_ggh4x = TRUE)
ggsave(plot = p1, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图sample_celltype.png', height = 7,width = 15,bg = 'white', dpi = 300,limitsize = FALSE)
ggsave(plot = p1, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图sample_celltype.pdf', height = 7,width = 15,bg = 'white', dpi = 300,limitsize = FALSE)

### Group4 3个分组
prop <- barplot_fraction(srt = scRNA,
                 celltype_col = "celltype",
                 group_col = "Group4",
                 sample_col = "sample")
if(is.numeric(prop$Freq)) prop$Freq <- 100 * prop$Freq else prop$Freq <- as.numeric(prop$Freq) * 100               
head(prop)
# 保存文件
write.table(prop, file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图sample_celltype_Group4.txt", sep = "\t", quote = F, row.names = F)

p1 <- plot_alluvial_barplot(data = prop, 
                      x = 'sample', 
                      y = 'Freq', 
                      fill_var = 'celltype',
                      fill_colors = RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                      strip_fill = c('#dfb67d','#cb8eab','#869746','#4f8799'),                                      
                      stratum_var = 'celltype', 
                      alluvium_var  ='celltype',
                      split.by = "Group4",
                      use_ggh4x = TRUE)
ggsave(plot = p1, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图sample_celltype_Group4.png', height = 7,width = 12,bg = 'white', dpi = 300,limitsize = FALSE)
ggsave(plot = p1, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图sample_celltype_Group4.pdf', height = 7,width = 12,bg = 'white', dpi = 300,limitsize = FALSE)
#-------------------------------------------------------------------------------------------------------------
### 分组与celltype
### Group1 4个分组
scRNA$new_group <- scRNA$Group1
prop <- barplot_fraction(srt = scRNA,
                 celltype_col = "celltype",
                 group_col = "Group1",
                 sample_col = "new_group")
if(is.numeric(prop$Freq)) prop$Freq <- 100 * prop$Freq else prop$Freq <- as.numeric(prop$Freq) * 100               
head(prop)
# 保存文件
write.table(prop, file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图Group1_celltype.txt", sep = "\t", quote = F, row.names = F)
p1 <- plot_alluvial_barplot(data = prop, 
                      x = 'Group1', 
                      y = 'Freq', 
                      fill_var = 'celltype',
                      fill_colors = RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                      strip_fill = c('#dfb67d','#cb8eab','#869746','#4f8799'),                                      
                      stratum_var = 'celltype', 
                      alluvium_var  ='celltype',
                      split.by = "Group1",
                      use_ggh4x = FALSE)
ggsave(plot = p1, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图Group1_celltype.png', height = 7,width = 6,bg = 'white', dpi = 300,limitsize = FALSE)
ggsave(plot = p1, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图Group1_celltype.pdf', height = 7,width = 6,bg = 'white', dpi = 300,limitsize = FALSE)               

#-------------------------------------------------------------------------------------------------------------
### 分组与celltype
### Group4 3个分组
scRNA$new_group <- scRNA$Group4
prop <- barplot_fraction(srt = scRNA,
                 celltype_col = "celltype",
                 group_col = "Group4",
                 sample_col = "new_group")
if(is.numeric(prop$Freq)) prop$Freq <- 100 * prop$Freq else prop$Freq <- as.numeric(prop$Freq) * 100               
head(prop)
# 保存文件
write.table(prop, file = "/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图Group4_celltype.txt", sep = "\t", quote = F, row.names = F)
p1 <- plot_alluvial_barplot(data = prop, 
                      x = 'Group4', 
                      y = 'Freq', 
                      fill_var = 'celltype',
                      fill_colors = RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                      strip_fill = c('#dfb67d','#cb8eab','#869746','#4f8799'),                                      
                      stratum_var = 'celltype', 
                      alluvium_var  ='celltype',
                      split.by = "Group1",
                      use_ggh4x = FALSE)
ggsave(plot = p1, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图Group4_celltype.png', height = 7,width = 5,bg = 'white', dpi = 300,limitsize = FALSE)
ggsave(plot = p1, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/堆叠柱形图Group4_celltype.pdf', height = 7,width = 5,bg = 'white', dpi = 300,limitsize = FALSE)  


#-------------------------------------------------------------------------------------------------------------
#                            查看不同基因的表达情况
#-------------------------------------------------------------------------------------------------------------
#读入不同分组过滤后的数据

#细胞增殖和生长调控基因
library(pheatmap)
library(ComplexHeatmap)

# 获得表达值
features <- c('AKT1', 'AKT2','CDK6', 'CCND1', 'CDK6', 'CDK6', 'CCNE1', 'CCNE2', 'MYC')
scRNA$new_group <- paste0(scRNA$celltype,'--',scRNA$Group1)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_39 <- t(df %>% dplyr::select(contains("37-39")))
rownames(df_37_39) <- stringr::str_split_fixed(rownames(df_37_39),pattern = "--",n = 2)[,1]
df_39_44 <- t(df %>% dplyr::select(contains("39-44")))
rownames(df_39_44) <- stringr::str_split_fixed(rownames(df_39_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 5, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_39,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_39',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_39_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '39_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p4 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞增殖和生长调控基因表达.png',height = 3,width = 15,units = 'in', res=300)
p1 + p2 + p3 + p4 
dev.off()

#****************************************************************************
scRNA$new_group <- paste0(scRNA$celltype,'--',scRNA$Group4)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_44 <- t(df %>% dplyr::select(contains("37-44")))
rownames(df_37_44) <- stringr::str_split_fixed(rownames(df_37_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 5, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞增殖和生长调控基因表达_Group4.png',height = 3,width = 15,units = 'in', res=300)
p1 + p2 + p3 
dev.off()

#****************************************************************************
scRNA$celltype_new <- ifelse(scRNA$celltype %in% c('CD56+CD16+NK','CD56+CD16-NK'),'NK cells','T cells')
scRNA$new_group <- paste0(scRNA$celltype_new,'--',scRNA$Group4)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_44 <- t(df %>% dplyr::select(contains("37-44")))
rownames(df_37_44) <- stringr::str_split_fixed(rownames(df_37_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 5, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞增殖和生长调控基因表达_Group4.png',height = 2,width = 12,units = 'in', res=300)
p1 + p2 + p3 
dev.off()



#-------------------------------------------------------------------------------------------------------------
#细胞周期基因
features <- c('CDC25A', 'CDK1','CDK2', 'CDK4', 'CDK6', 'FOXO1', 'CDKN2A', 'CDKN2B')
scRNA$new_group <- paste0(scRNA$celltype,'--',scRNA$Group1)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_39 <- t(df %>% dplyr::select(contains("37-39")))
rownames(df_37_39) <- stringr::str_split_fixed(rownames(df_37_39),pattern = "--",n = 2)[,1]
df_39_44 <- t(df %>% dplyr::select(contains("39-44")))
rownames(df_39_44) <- stringr::str_split_fixed(rownames(df_39_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 5, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_39,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_39',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_39_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '39_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p4 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞周期基因.png',height = 4,width = 18,units = 'in', res=300)
p1 + p2 + p3 + p4 
dev.off()

#****************************************************************************
scRNA$new_group <- paste0(scRNA$celltype,'--',scRNA$Group4)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_44 <- t(df %>% dplyr::select(contains("37-44")))
rownames(df_37_44) <- stringr::str_split_fixed(rownames(df_37_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 5, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞周期基因_Group4.png',height = 3,width = 12,units = 'in', res=300)
p1 + p2 + p3 
dev.off()

#****************************************************************************
scRNA$celltype_new <- ifelse(scRNA$celltype %in% c('CD56+CD16+NK','CD56+CD16-NK'),'NK cells','T cells')
scRNA$new_group <- paste0(scRNA$celltype_new,'--',scRNA$Group4)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_44 <- t(df %>% dplyr::select(contains("37-44")))
rownames(df_37_44) <- stringr::str_split_fixed(rownames(df_37_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 5, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞周期基因_大群.png',height = 2,width = 12,units = 'in', res=300)
p1 + p2 + p3 
dev.off()

#****************************************************************************
scRNA$new_group <- paste0(scRNA$sample,'--',scRNA$celltype,'--',scRNA$Group4)
df <- t(as.data.frame(AverageExpression(scRNA, features = 'CDK6',slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA))
df <- as.data.frame(df)
df$sample <- stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,1]
df$celltype <- stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,2]
df$group = stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,3]
colnames(df)[1] <- 'Expression'

p <- ggplot(df,mapping = aes(x=group,y=Expression,fill=group))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.2,size=1) +
    facet_wrap(vars(celltype),ncol=5)+
    labs(x = '',y='Sample Average Expression') + 
    ggtitle("CDK6")+
    theme(plot.title = element_text(hjust = 0.5,size=18,color = 'black', face='bold'),
          strip.text = element_text(size=14,color='black'),
          axis.text.x = element_text(size = 12,angle = 90,hjust = 0.5,vjust = 0.5))

ggsave(plot = p, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/基因表达/001_celltype_Group1_CDK6.png', height = 4,width = 10,bg = 'white', dpi = 200,limitsize = FALSE)


#*****************************************************************************
scRNA$new_group <- paste0(scRNA$sample,'--',scRNA$celltype,'--',scRNA$Group1)
df <- t(as.data.frame(AverageExpression(scRNA, features = 'CDK6',slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA))
df <- as.data.frame(df)
df$sample <- stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,1]
df$celltype <- stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,2]
df$group = stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,3]
colnames(df)[1] <- 'Expression'

p <- ggplot(df,mapping = aes(x=group,y=Expression,fill=group))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.2,size=1) +
    facet_wrap(vars(celltype),ncol=5)+
    labs(x = '',y='Sample Average Expression') + 
    ggtitle("CDK6")+
    theme(plot.title = element_text(hjust = 0.5,size=18,color = 'black', face='bold'),
          strip.text = element_text(size=14,color='black'),
          axis.text.x = element_text(size = 12,angle = 90,hjust = 0.5,vjust = 0.5))

ggsave(plot = p, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/基因表达/celltype_Group2_CDK6.png', height = 4,width = 10,bg = 'white', dpi = 200,limitsize = FALSE)



#-------------------------------------------------------------------------------------------------------------
#                            绘制FOXO1的基因表达箱线图
#-------------------------------------------------------------------------------------------------------------
scRNA$new_group <- paste0(scRNA$sample,'--',scRNA$celltype,'--',scRNA$Group4)
df <- t(as.data.frame(AverageExpression(scRNA, features = 'FOXO1',slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA))
df <- as.data.frame(df)
df$sample <- stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,1]
df$celltype <- stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,2]
df$group = stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,3]
colnames(df)[1] <- 'Expression'

p <- ggplot(df,mapping = aes(x=group,y=Expression,fill=group))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.2,size=1) +
    facet_wrap(vars(celltype),ncol=5)+
    labs(x = '',y='Sample Average Expression') + 
    ggtitle("FOXO1")+
    theme(plot.title = element_text(hjust = 0.5,size=18,color = 'black', face='bold'),
          strip.text = element_text(size=14,color='black'),
          axis.text.x = element_text(size = 12,angle = 90,hjust = 0.5,vjust = 0.5))

ggsave(plot = p, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/基因表达/celltype_Group1_FOXO1.png', height = 4,width = 10,bg = 'white', dpi = 200,limitsize = FALSE)


#*****************************************************************************
scRNA$new_group <- paste0(scRNA$sample,'--',scRNA$celltype,'--',scRNA$Group1)
df <- t(as.data.frame(AverageExpression(scRNA, features = 'FOXO1',slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA))
df <- as.data.frame(df)
df$sample <- stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,1]
df$celltype <- stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,2]
df$group = stringr::str_split_fixed(rownames(df),pattern = "--",n = 3)[,3]
colnames(df)[1] <- 'Expression'

p <- ggplot(df,mapping = aes(x=group,y=Expression,fill=group))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.2,size=1) +
    facet_wrap(vars(celltype),ncol=5)+
    labs(x = '',y='Sample Average Expression') + 
    ggtitle("FOXO1")+
    theme(plot.title = element_text(hjust = 0.5,size=18,color = 'black', face='bold'),
          strip.text = element_text(size=14,color='black'),
          axis.text.x = element_text(size = 12,angle = 90,hjust = 0.5,vjust = 0.5))

ggsave(plot = p, file = '/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/基因表达/celltype_Group2_FOXO1.png', height = 4,width = 10,bg = 'white', dpi = 200,limitsize = FALSE)




#-------------------------------------------------------------------------------------------------------------
# 信号传导和代谢调节基因
features <- c('PIK3CA', 'PIK3CB', 'PIK3CD', 'FOXO1', 'PIK3R2', 'PIK3R3', 'KRAS', 'NRAS', 'HRAS', 'FOXO1', 'FOXO3', 'MTOR','MAP2K1', 'MAP2K2', 'MAP2K3', 'MAP2K6', 'MAPK1', 'MAPK11',
              'MAPK12', 'MAPK13', 'MAPK14', 'MAPK3','MAPKAPK2', 'RAF1', 'RASSF5', 'RRAS', 'RRAS2', 'RHEB', 'MRAS', 'HUS1')
scRNA$new_group <- paste0(scRNA$celltype,'--',scRNA$Group1)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_39 <- t(df %>% dplyr::select(contains("37-39")))
rownames(df_37_39) <- stringr::str_split_fixed(rownames(df_37_39),pattern = "--",n = 2)[,1]
df_39_44 <- t(df %>% dplyr::select(contains("39-44")))
rownames(df_39_44) <- stringr::str_split_fixed(rownames(df_39_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 5, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_39,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_39',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_39_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '39_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p4 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/基因表达/信号传导和代谢调节基因.png',height = 3,width = 35,units = 'in', res=300)
p1 + p2 + p3 + p4 
dev.off()

#******************************************************
# 基因集打分
features_list <- list(single = features)
scRNA <- Seurat::AddModuleScore(scRNA, features = features_list, name = 'single')
df <- scRNA@meta.data %>% group_by(sample,Group1) %>% dplyr::summarise(mean = mean(single1))
p1 <- ggplot(df,mapping = aes(x=Group1,y=mean,fill=Group1))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.2)+
    labs(x = '',y='Sample Average Expression') + 
    # ggtitle("FOXO1")+
    theme(plot.title = element_text(hjust = 0.5,size=18,color = 'black', face='bold'),
          strip.text = element_text(size=14,color='black'),
          axis.text.x = element_text(size = 12,angle = 90,hjust = 0.5,vjust = 0.5))

df2 <- scRNA@meta.data %>% group_by(sample,Group4) %>% dplyr::summarise(mean = mean(single1))
p2 <- ggplot(df2,mapping = aes(x=Group4,y=mean,fill=Group4))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.2)+
    labs(x = '',y='Sample Average Expression') + 
    # ggtitle("FOXO1")+
    theme(plot.title = element_text(hjust = 0.5,size=18,color = 'black', face='bold'),
          strip.text = element_text(size=14,color='black'),
          axis.text.x = element_text(size = 12,angle = 90,hjust = 0.5,vjust = 0.5))

p1+p2+ plot_layout(ncol = 2)+plot_annotation(tag_levels = 'A')


#-------------------------------------------------------------------------------------------------------------
# 抑制细胞衰老的基因
features <- c('PTEN', 'TP53','RB1', 'RBL1', 'RBL2')
scRNA$new_group <- paste0(scRNA$celltype,'--',scRNA$Group1)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_39 <- t(df %>% dplyr::select(contains("37-39")))
rownames(df_37_39) <- stringr::str_split_fixed(rownames(df_37_39),pattern = "--",n = 2)[,1]
df_39_44 <- t(df %>% dplyr::select(contains("39-44")))
rownames(df_39_44) <- stringr::str_split_fixed(rownames(df_39_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 3, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_39,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_39',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_39_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '39_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p4 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/基因表达/抑制细胞衰老的基因.png',height = 3,width = 10,units = 'in', res=300)
p1 + p2 + p3 + p4 
dev.off()











#-------------------------------------------------------------------------------------------------------------
# 细胞衰老诱导基因
features <-c('GADD45A', 'GADD45B', 'GADD45G', 'SERPINE1', 'TGFB1', 'TGFB2', 'TGFB3')
scRNA$new_group <- paste0(scRNA$celltype,'--',scRNA$Group1)
df <- as.data.frame(AverageExpression(scRNA, features = features,slot = 'data',group.by = 'new_group',return.seurat = FALSE)$RNA)
df_18_35 <- t(df %>% dplyr::select(contains("18-35")))
rownames(df_18_35) <- stringr::str_split_fixed(rownames(df_18_35),pattern = "--",n = 2)[,1]
df_37_39 <- t(df %>% dplyr::select(contains("37-39")))
rownames(df_37_39) <- stringr::str_split_fixed(rownames(df_37_39),pattern = "--",n = 2)[,1]
df_39_44 <- t(df %>% dplyr::select(contains("39-44")))
rownames(df_39_44) <- stringr::str_split_fixed(rownames(df_39_44),pattern = "--",n = 2)[,1]
df_47_49 <- t(df %>% dplyr::select(contains("47-49")))
rownames(df_47_49) <- stringr::str_split_fixed(rownames(df_47_49),pattern = "--",n = 2)[,1]

# 设置颜色标度
library(circlize)
# colorRampPalette(rev(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50)) ##蓝到红
colors<-colorRampPalette(rev(c("#B2182B","#EF8A62","white","#67A9CF","#2166AC")))(50)
values <- seq(0, 5, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)
p1 <- Heatmap(df_18_35,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '18_35',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p2 <- Heatmap(df_37_39,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '37_39',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p3 <- Heatmap(df_39_44,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns  = FALSE,column_title = '39_44',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))
p4 <- Heatmap(df_47_49,col = col_fun, name = 'Expression',cluster_rows = FALSE, cluster_columns =  FALSE,column_title = '47_49',column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12))

png('/root/wangje/Project/OvaryFigs/02_NKandTcells/UMAP/大群/细胞衰老诱导基因.png',height = 5,width = 15,units = 'in', res=300)
p1 + p2 + p3 + p4 
dev.off()

#-------------------------------------------------------------------------------------------------------------
#                            富集分析
#-------------------------------------------------------------------------------------------------------------


