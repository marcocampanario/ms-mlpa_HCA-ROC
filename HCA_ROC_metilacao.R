### SCRIPT para clusterizacao de genes com base em seu padrao de metilacao

### Bioinformata: Marco Antonio Campanario (#macs)
### Contato: macscampanario@gmail.com
### Ultima atualizacao: 22.04.2025

library(grid)
library(openxlsx)
library(pheatmap)  # para visualizacao
library(cluster)   # para analise
library(factoextra)  # visualizacao elegante
library(tidyverse)
library(ggplot2)
library(plotROC)
library(pROC)

data <- readxl::read_excel("ms-mlpa_tumor_adjacent.xlsx")

# criar uma TAG para tecido TUMORAL (tum) e ADJACENTE (adj)

data_tecido <- data %>%
  mutate(Tecido = c(rep("tum", 81), rep("adj", 15))) %>%
  select(Tecido, everything())

# deixar somente colunas com ID, Tecido e Porcentagens de Metilacao por gene

data_tecido_pct <- data_tecido %>%
  select(names(data_tecido)[c(1, 2, 23:83)]) %>%
  select(-matches("_CAT")) %>%
  select(-c(31:34))

# TUMOR X ADJACENTE ------------------------------------------------------------

## Clusterizacao por porcentagem de metilacao TUMOR X ADJACENTE

# retirar os genes fora do kit (tag 999)

data_tecido_pct <- data_tecido_pct %>%
  select(where(~ !any(. == 999))) # 2 genes sairam (PTEN e KLLN)

# filtrar apenas os pacientes que tem tecido adjacente e tumoral
data_pareado <- data_tecido_pct %>%
  filter(ID %in% 2:18)

# separar os dados de tumor e normal
tumor <- data_pareado %>% filter(Tecido == "tum") %>% arrange(ID)
adjacente <- data_pareado %>% filter(Tecido == "adj") %>% arrange(ID)

# verificar se estao pareados
identical(tumor$ID, adjacente$ID)

# calcular a diferença de metilação: Tumor - Normal
diferencas <- tumor[, 3:ncol(tumor)] - adjacente[, 3:ncol(adjacente)]
rownames(diferencas) <- paste0("Paciente_", tumor$ID)

# verificar se tem valores NA ou NaN nos dados que vao para clusterização
any(is.na(diferencas))
any(is.nan(as.matrix(diferencas)))
any(is.infinite(as.matrix(diferencas)))

# verificar se tem colunas constantes e retirar
apply(diferencas, 2, function(x) length(unique(x)))
diferencas <- diferencas[, apply(diferencas, 2, function(x) length(unique(x)) > 1)]
# 3 genes sairam (KLLN_A, VHL e APC)

# transpor o dataframe (para clusterizar genes e nao pacientes)
diferencas_t <- t(diferencas)

# calcular a matriz de distâncias
dist_mat <- dist(diferencas_t, method = "euclidean")

# cluster hierárquico
hc <- hclust(dist_mat, method = "ward.D2")

# dendrograma
plot(hc, main = "clusterizacao_hierarquica_genes_TESTE", xlab = "", sub = "")

# heatmap
pheatmap(diferencas_t,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         main = "Padrão de metilação (Tumor - Normal)")

# pca
pca <- prcomp(diferencas_t, scale. = TRUE)
fviz_pca_biplot(pca, repel = TRUE)

# SOMENTE TUMOR ----------------------------------------------------------------

## Clusterizacao por porcentagem de metilacao SOMENTE TUMORES

# selecionar as amostras tumorais
tumores <- data_tecido_pct %>%
  filter(Tecido == "tum") %>%
  arrange(ID)

# substituir "999" por NA (valores ausentes)
tumores[tumores == 999] <- NA

# remover colunas ID e Tecido
tumores_clean <- tumores %>%
  select(-c(ID, Tecido))

# transpor para clusterização por genes
tumores_t <- t(tumores_clean)

# guardar matriz com NA para o heatmap
tumores_t_heatmap <- tumores_t

# substituir NA por média da linha (gene), só para o dist/hclust
tumores_t_nona <- t(apply(tumores_t, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))

# calculo de matriz de distancia com dados sem NA

# distancia geometrica (euclidean)
dist_mat <- dist(tumores_t_nona, method = "euclidean")
hc_mat <- hclust(dist_mat, method = "ward.D2")
plot(hc_mat, main = "Dendrograma da clusterização dos genes por distância geométrica")

# heatmap usando a matriz original com NA, mas cluster definido pelos dados sem NA (euclidiana)
pheatmap(tumores_t_heatmap,
         cluster_rows = hc_mat,
         cluster_cols = TRUE,
         main = "Clusteriza??o hier?rquica do perfil de metila??o (%) em tumores colorretais",
         scale = "none",
         color = colorRampPalette(c("black", "red", "green"))(100),
         na_col = "grey",
         legend = TRUE,
         breaks = seq(0, 100, length.out = 101),
         legend_breaks = c(0, 25, 50, 75, 100),
         legend_labels = c("0%", "25%", "50%", "75%", "100%"),
         border_color = NA)


# CURVA ROC --------------------------------------------------------------------

# retirar a coluna de ID das amostras
roc <- data_tecido_pct %>%
  select(-ID)

# trocar o código tum/adj para 1/0
roc$Tecido <- ifelse(roc$Tecido == "tum", 1,
                        ifelse(roc$Tecido == "adj", 0, NA))

# trocar o nome da coluna Tecido para D (padrão pacote pROC) 
names(roc)[names(roc) == "Tecido"] <- "D"

# cálculo ROC e plotagem da curva

## ESR1
# cálculo ROC + especificidade/sensibilidade + melhor ponto de corte
roc_obj_ESR1 <- roc(response = roc$D, predictor = roc$ESR1)
coords_all_ESR1 <- coords(roc_obj_ESR1, x = "all", ret = c("threshold", "sensitivity", "specificity",
                                                             "accuracy",
                                                             "ppv",
                                                             "npv",
                                                             "fpr",
                                                             "fnr"),
                           transpose = F)
coords_best_ESR1 <- coords(roc_obj_ESR1, x = "best", best.method = "youden", ret = c("threshold",
                                                                                     "sensitivity",
                                                                                     "specificity"))
# plotagem
ESR1 <- ggplot(roc, aes(d = D, m = ESR1)) +
  geom_roc(n.cuts = 0) +
  geom_rocci(ci.at = as.numeric(quantile(roc$ESR1, c(.1, .4, .5, .6, .9), na.rm = TRUE))) +
  style_roc()
direct_label(ESR1, labels = "ESR1", nudge_y = -.2)
calc_auc(ESR1)

## CDH13
# cálculo ROC + especificidade/sensibilidade + melhor ponto de corte
roc_obj_CDH13 <- roc(response = roc$D, predictor = roc$CDH13)
coords_all_CDH13 <- coords(roc_obj_CDH13, x = "all", ret = c("threshold", "sensitivity", "specificity",
                                                 "accuracy",
                                                 "ppv",
                                                 "npv",
                                                 "fpr",
                                                 "fnr"),
                     transpose = F)
coords_best_CDH13 <- coords(roc_obj_CDH13, x = "best", best.method = "youden", ret = c("threshold",
                                                                                       "sensitivity",
                                                                                       "specificity"))
# plotagem
CDH13 <- ggplot(roc, aes(d = D, m = CDH13)) +
  geom_roc(n.cuts = 0) +
  geom_rocci(ci.at = as.numeric(quantile(roc$CDH13, c(.1, .4, .5, .6, .9), na.rm = TRUE))) +
  style_roc()
direct_label(CDH13, labels = "CDH13", nudge_y = -.2)

## RARB
# cálculo ROC + especificidade/sensibilidade + melhor ponto de corte
roc_obj_RARB <- roc(response = roc$D, predictor = roc$RARB)
coords_all_RARB <- coords(roc_obj_RARB, x = "all", ret = c("threshold", "sensitivity", "specificity",
                                                           "accuracy",
                                                           "ppv",
                                                           "npv",
                                                           "fpr",
                                                           "fnr"),
                          transpose = F)
coords_best_RARB <- coords(roc_obj_RARB, x = "best", best.method = "youden", ret = c("threshold",
                                                                                     "sensitivity",
                                                                                     "specificity"))
# plotagem
RARB <- ggplot(roc, aes(d = D, m = RARB)) +
  geom_roc(n.cuts = 0) +
  geom_rocci(ci.at = as.numeric(quantile(roc$RARB, c(.1, .4, .5, .6, .9), na.rm = TRUE))) +
  style_roc()
direct_label(RARB, labels = "RARB", nudge_y = -.2)
calc_auc(RARB)

## CHFR
# cálculo ROC + especificidade/sensibilidade + melhor ponto de corte
roc_obj_CHFR <- roc(response = roc$D, predictor = roc$CHFR)
coords_all_CHFR <- coords(roc_obj_CHFR, x = "all", ret = c("threshold", "sensitivity", "specificity",
                                                           "accuracy",
                                                           "ppv",
                                                           "npv",
                                                           "fpr",
                                                           "fnr"),
                          transpose = F)
coords_best_CHFR <- coords(roc_obj_CHFR, x = "best", best.method = "youden", ret = c("threshold",
                                                                                     "sensitivity",
                                                                                     "specificity"))
# plotagem
CHFR <- ggplot(roc, aes(d = D, m = CHFR)) +
  geom_roc(n.cuts = 0) +
  geom_rocci(ci.at = as.numeric(quantile(roc$CHFR, c(.1, .4, .5, .6, .9), na.rm = TRUE))) +
  style_roc()
direct_label(CHFR, labels = "CHFR", nudge_y = -.2)

## p-valores (Teste U de Mann-Whitney / Wilcoxon rank-sum test)
## Nao recomendado para dados com muitos valores identicos (0)

genes <- c("ESR1", "CDH13", "RARB", "CHFR")
resultados <- data.frame(Gene = character(), AUC = numeric(), P_Value = numeric())
for (g in genes) {
  predictor <- roc[[g]]  # met %
  group <- roc$D  # 0 = normal, 1 = tumor
  
  # ROC/AUC
  r <- roc(group, predictor)
  auc_value <- auc(r)
  
  # Wilcoxon rank-sum test
  wilcox_res <- wilcox.test(predictor ~ group)
  p_val <- wilcox_res$p.value

  resultados <- rbind(resultados, data.frame(genes = g, AUC = as.numeric(auc_value), P_Value = p_val))
}
print(resultados)

## Alternativa: Teste de Permutacao

# AUC
calc_auc <- function(true_labels, predictions) {
  r <- roc(true_labels, predictions)
  return(auc(r))
}

# Permutacoes
n_permutations <- 100000

results <- data.frame(Gene = character(), AUC = numeric(), P_Value = numeric())

# Genes
genes <- c("ESR1", "CDH13", "RARB", "CHFR")

# Loop
for (g in genes) {
  
  # Met %
  predictor <- roc[[g]]
  
  # AUC Observada
  observed_auc <- calc_auc(roc$D, predictor)
  
  # Permutacoes
  permuted_aucs <- numeric(n_permutations)
  set.seed(42)  # Reprodutibilidade
  for (i in 1:n_permutations) {
    permuted_labels <- sample(roc$D)  # Aleatorizacao
    permuted_aucs[i] <- calc_auc(permuted_labels, predictor)  # AUC permutada
  }
  
  # P-valor: Proporcao de AUCs permutadas maiores ou iguais a AUC observada
  p_value <- mean(permuted_aucs >= observed_auc)
  
  results <- rbind(results, data.frame(genes = g, AUC = observed_auc, P_Value = p_value))
}

#
print(results)
                                 
## Visualizacao de todos os 4

# criar objeto novo para trabalhar a visualizacao conjunta
genes_roc <- roc

# fazer o melt da tabela com a metilacao de cada um dos 4 genes
genes_roc <- melt_roc(genes_roc, "D", c("ESR1", "CDH13", "RARB", "CHFR"))

# criar o plot da visualizacao conjunta
genes_plot <- ggplot(genes_roc, aes(d = D.D, m = M)) + geom_roc() + facet_wrap(~ name) + style_roc() +
  theme(
      strip.text = element_text(size = 16, face = "bold"),     
      axis.title.x = element_text(size = 12, face = "bold"),   
      axis.title.y = element_text(size = 12, face = "bold")
    )

# calcular a area sob curva/melhor limiar de corte de cada gene e
# plotar na visualizacao conjunta
auc <- genes_roc %>%
  group_by(name) %>%
  group_modify(~ {
    roc_plot <- ggplot(.x, aes(d = D.D, m = M)) + geom_roc(n.cuts = 0)
    calc_auc(roc_plot)
  }) %>%
  ungroup()
auc <- auc %>%
  mutate(label = paste0("AUC = ", round(AUC, 3)))

# tabela melhor threshold dos 4
best_threshold <- data.frame(name = c("ESR1", "CDH13", "RARB", "CHFR"),
                             best_threshold = c(coords_best_ESR1$threshold,
                                              coords_best_CDH13$threshold,
                                              coords_best_RARB$threshold,
                                              coords_best_CHFR$threshold))
auc <- inner_join(auc, best_threshold, by = "name")
  
# plotagem da visualizacao conjunta 
genes_plot +
geom_text(data = auc,
            aes(x = 0.75, y = 0.25, label = paste0(label, "\n", "Best threshold = ", best_threshold)),
            inherit.aes = FALSE,
            size = 4)

#### ESTADO EMPIRICO DE METILACAO ----------------------------------------------

data_empirico <- data

limiares <- list(
  "CDH13" = coords_best_CDH13$threshold,
  "CHFR" = coords_best_CHFR$threshold,
  "ESR1" = coords_best_ESR1$threshold,
  "RARB" = coords_best_RARB$threshold
)

for (gene in names(limiares)) {
  # acha a coluna com os dados de metilacao
  idx <- grep(paste0("^", gene, "$"), names(data_empirico))
  
  if (length(idx) == 1) {
    col_met <- data_empirico[[idx]]
    limiar <- limiares[[gene]]
    estado_empirico <- ifelse(col_met > limiar, 2, 1)
    
    # adiciona a nova coluna ao lado direito da coluna original
    nova_coluna <- data.frame(estado_empirico)
    colnames(nova_coluna) <- paste0(gene, "_empirico")
    data_empirico <- cbind(
      data_empirico[ , 1:idx],
      nova_coluna,
      if (idx < ncol(data_empirico)) data_empirico[ , (idx + 1):ncol(data_empirico)] else NULL
    )
  } else {
    warning(paste("Gene", gene, "não encontrado ou encontrado mais de uma vez."))
  }
}

# FIGURAS E TABELAS ------------------------------------------------------------

# DENDROGRAMA
tiff("dendrograma_genes_en.tiff", width = 10, height = 8, units = "in", res = 300)
plot(hc_mat, main = "Dendrogram of gene clustering by geometric distance", sub = "", xlab = "", ylab = "Distance")
dev.off()

tiff("dendrograma_genes_ptbr.tiff", width = 10, height = 8, units = "in", res = 300)
plot(hc_mat, main = "Dendrograma da clusterização dos genes por distância geométrica", sub = "", xlab = "", ylab = "Distância")
dev.off()

# HEATMAP
tiff("heatmap_tumores_ptbr.tiff", width = 10, height = 10, units = "in", res = 300)
pheatmap(tumores_t_heatmap,
         cluster_rows = hc_mat,
         cluster_cols = TRUE,
         main = "Clusterização hierárquica do perfil de metilação (%) em tumores colorretais",
         scale = "none",
         color = colorRampPalette(c("black", "red", "green"))(100),
         na_col = "grey",
         legend = TRUE,
         breaks = seq(0, 100, length.out = 101),
         legend_breaks = c(0, 25, 50, 75, 100),
         legend_labels = c("0%", "25%", "50%", "75%", "100%"),
         border_color = NA)
dev.off()

tiff("heatmap_tumores_en.tiff", width = 10, height = 10, units = "in", res = 300)
pheatmap(tumores_t_heatmap,
         cluster_rows = hc_mat,
         cluster_cols = TRUE,
         main = "Hierarchical clustering of methylation profile (%) in colorectal tumors",
         scale = "none",
         color = colorRampPalette(c("black", "red", "green"))(100),
         na_col = "grey",
         legend = TRUE,
         breaks = seq(0, 100, length.out = 101),
         legend_breaks = c(0, 25, 50, 75, 100),
         legend_labels = c("0%", "25%", "50%", "75%", "100%"),
         border_color = NA)
dev.off()

# TABELA SENSIBILIDADE/ESPECIFICIDADE ESR1
coords_all_ESR1
# TABELA SENSIBILIDADE/ESPECIFICIDADE CDH13
coords_all_CDH13
# TABELA SENSIBILIDADE/ESPECIFICIDADE RARB
coords_all_RARB
# TABELA SENSIBILIDADE/ESPECIFICIDADE CHFR
coords_all_CHFR
# TABELA AREA SOB CURVA + MELHOR LIMIAR DE CORTE
AUC <- auc %>%
  select(c("name", "AUC", "best_threshold"))
names(AUC) <- c("genes", "auc", "best_threshold")
# TABELA DADOS EMPIRICOS
data_empirico

# TABELA AGREGADA
wb <- createWorkbook()
addWorksheet(wb, "ESR1")
writeData(wb, "ESR1", coords_all_ESR1)
addWorksheet(wb, "CDH13")
writeData(wb, "CDH13", coords_all_CDH13)
addWorksheet(wb, "RARB")
writeData(wb, "RARB", coords_all_RARB)
addWorksheet(wb, "CHFR")
writeData(wb, "CHFR", coords_all_CHFR)
addWorksheet(wb, "AUC_best_threshold")
writeData(wb, "AUC_best_threshold", data_empirico)
addWorksheet(wb, "data_empirico")
writeData(wb, "data_empirico", data_empirico)
saveWorkbook(wb, "tabelas_roc.xlsx", overwrite = TRUE)


# CURVAS ROC - VISUALIZACAO CONJUNTA
tiff("curvas_roc.tiff", width = 10, height = 10, units = "in", res = 300)
genes_plot <- ggplot(genes_roc, aes(d = D.D, m = M)) +
  geom_roc() +
  facet_wrap(~ name) +
  style_roc() +  
  labs(
    x = "False Positive Fraction",
    y = "True Positive Fraction"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.07, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.05))) +
  coord_cartesian(clip = "off") +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 30, 30)
  )
genes_plot +
  geom_text(data = auc,
            aes(x = 0.75, y = 0.25, label = paste0(label, "\n", "Best threshold = ", best_threshold)),
            inherit.aes = FALSE,
            size = 4)
dev.off()





