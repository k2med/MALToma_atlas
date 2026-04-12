library(msigdbr)
library(dplyr)

genesets_HALLMARK <- msigdbr(species = "Homo sapiens"
                             , category = "H") %>%
  subset(select = c("gs_name","gene_symbol")) %>%
  as.data.frame() %>% unique()

genesets_GO_BP <- msigdbr(species = "Homo sapiens"
                          , category = "C5", subcategory = 'GO:BP') %>%
  subset(select = c("gs_name","gene_symbol")) %>%
  as.data.frame() %>% unique()

genesets_GO_MF <- msigdbr(species = "Homo sapiens"
                          , category = "C5", subcategory = 'GO:MF') %>%
  subset(select = c("gs_name","gene_symbol")) %>%
  as.data.frame() %>% unique()

genesets_GO_CC <- msigdbr(species = "Homo sapiens"
                          , category = "C5", subcategory = 'GO:CC') %>%
  subset(select = c("gs_name","gene_symbol")) %>%
  as.data.frame() %>% unique()

genesets_KEGG <- msigdbr(species = "Homo sapiens"
                         , category = "C2", subcategory = "KEGG") %>%
  subset(select = c("gs_name","gene_symbol")) %>%
  as.data.frame() %>% unique()

genesets_REACTOME <- msigdbr(species = "Homo sapiens"
                             , category = "C2", subcategory = 'CP:REACTOME') %>%
  subset(select = c("gs_name","gene_symbol")) %>%
  as.data.frame() %>% unique()

genesets_merge <- rbind(genesets_HALLMARK, genesets_GO_BP, genesets_GO_MF,
                        genesets_GO_CC, genesets_KEGG, genesets_REACTOME)

colnames(genesets_merge) <- c('pathID','geneID')

head(genesets_merge)

format_name <- function(x) {
  prefix <- sub("^(.*?)_.*", "\\1", x)       # 获取前缀
  body <- sub("^[^_]+_", "", x)              # 移除前缀
  words <- tolower(unlist(strsplit(body, "_"))) # 转小写并拆分
  # 将首词首字母大写，不依赖 tools::toTitleCase
  words[1] <- paste0(toupper(substr(words[1], 1, 1)), substr(words[1], 2, nchar(words[1])))
  formatted <- paste(words, collapse = " ")  # 合并
  paste0(formatted, " [", prefix, "]")       # 添加后缀
}
sapply(unique(genesets_merge$pathID)[1:100],  format_name)

genesets_merge$pathID_rename <- sapply(genesets_merge$pathID, format_name)

saveRDS(genesets_merge, 'msigdbr_gene_sets.rds')
