# ============================================================
# 03 - ESCALA DE EFEITO
# Veado-catingueiro - curva BIMODAL (escolhida no script 02)
#
# Para cada variavel de paisagem, ajusta o modelo nas 6 escalas
# (800, 1240, 1680, 2120, 2560, 3000 m) de duas formas:
#   aditivo   : time + x  -> a variavel muda so o NIVEL de atividade
#   interacao : time * x  -> a variavel muda tambem o HORARIO
# A melhor escala e a de menor AIC. O modelo nulo (so time) entra em
# cada tabela, ajustado nos MESMOS dados.
#
# Variaveis (padronizadas, prefixo z_; sao do nivel de cluster):
#   A. native    - % vegetacao nativa
#   B. ndvi_mean - NDVI medio da vegetacao nativa
#   C. enn_mn    - distancia media entre fragmentos
#
# Entrada: data/processed/chm_veado.rds
# Saidas : output/tabelas/03_*.xlsx, output/figuras/03_escala_efeito.png
# ============================================================


# ------------------------------------------------------------
# 0. Pacotes, funcoes do CHM e dados
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(openxlsx)
library(here)
source(here("R", "chm_helpers.R"))

dir.create(here("output", "tabelas"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figuras"), recursive = TRUE, showWarnings = FALSE)

chm_veado <- readRDS(here("data", "processed", "chm_veado.rds"))

# Modelo nulo (mesmo do script 02): bimodal, cameras aninhadas em clusters
nulo <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
               data = chm_veado)

summary(nulo)
# ============================================================
# A. VEGETACAO NATIVA (native)
# ============================================================
native_800_aditivo  <- chmTMB(cbind(success, failure) ~ time + z_native_800  + (1 | cluster / locationName), data = chm_veado)
native_1240_aditivo <- chmTMB(cbind(success, failure) ~ time + z_native_1240 + (1 | cluster / locationName), data = chm_veado)
native_1680_aditivo <- chmTMB(cbind(success, failure) ~ time + z_native_1680 + (1 | cluster / locationName), data = chm_veado)
native_2120_aditivo <- chmTMB(cbind(success, failure) ~ time + z_native_2120 + (1 | cluster / locationName), data = chm_veado)
native_2560_aditivo <- chmTMB(cbind(success, failure) ~ time + z_native_2560 + (1 | cluster / locationName), data = chm_veado)
native_3000_aditivo <- chmTMB(cbind(success, failure) ~ time + z_native_3000 + (1 | cluster / locationName), data = chm_veado)

native_800_interacao  <- chmTMB(cbind(success, failure) ~ time * z_native_800  + (1 | cluster / locationName), data = chm_veado)
native_1240_interacao <- chmTMB(cbind(success, failure) ~ time * z_native_1240 + (1 | cluster / locationName), data = chm_veado)
native_1680_interacao <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName), data = chm_veado)
native_2120_interacao <- chmTMB(cbind(success, failure) ~ time * z_native_2120 + (1 | cluster / locationName), data = chm_veado)
native_2560_interacao <- chmTMB(cbind(success, failure) ~ time * z_native_2560 + (1 | cluster / locationName), data = chm_veado)
native_3000_interacao <- chmTMB(cbind(success, failure) ~ time * z_native_3000 + (1 | cluster / locationName), data = chm_veado)

cat("\n===== A. VEGETACAO NATIVA =====\n")
cat("Convergencia (TRUE = ok):\n")
print(c(nulo = nulo$sdr$pdHess,
        native_800_aditivo  = native_800_aditivo$sdr$pdHess,
        native_1240_aditivo = native_1240_aditivo$sdr$pdHess,
        native_1680_aditivo = native_1680_aditivo$sdr$pdHess,
        native_2120_aditivo = native_2120_aditivo$sdr$pdHess,
        native_2560_aditivo = native_2560_aditivo$sdr$pdHess,
        native_3000_aditivo = native_3000_aditivo$sdr$pdHess,
        native_800_interacao  = native_800_interacao$sdr$pdHess,
        native_1240_interacao = native_1240_interacao$sdr$pdHess,
        native_1680_interacao = native_1680_interacao$sdr$pdHess,
        native_2120_interacao = native_2120_interacao$sdr$pdHess,
        native_2560_interacao = native_2560_interacao$sdr$pdHess,
        native_3000_interacao = native_3000_interacao$sdr$pdHess))

tab_native <- AIC(nulo,
                  native_800_aditivo, native_1240_aditivo, native_1680_aditivo,
                  native_2120_aditivo, native_2560_aditivo, native_3000_aditivo,
                  native_800_interacao, native_1240_interacao, native_1680_interacao,
                  native_2120_interacao, native_2560_interacao, native_3000_interacao) %>%
  tibble::rownames_to_column("modelo") %>%
  mutate(variavel = "native",
         n_cameras = n_distinct(chm_veado$locationName))
print(tab_native %>% mutate(dAIC = round(AIC - min(AIC), 2)) %>% arrange(dAIC))


# ============================================================
# B. NDVI MEDIO (ndvi_mean)
# ============================================================
ndvi_800_aditivo  <- chmTMB(cbind(success, failure) ~ time + z_ndvi_mean_800  + (1 | cluster / locationName), data = chm_veado)
ndvi_1240_aditivo <- chmTMB(cbind(success, failure) ~ time + z_ndvi_mean_1240 + (1 | cluster / locationName), data = chm_veado)
ndvi_1680_aditivo <- chmTMB(cbind(success, failure) ~ time + z_ndvi_mean_1680 + (1 | cluster / locationName), data = chm_veado)
ndvi_2120_aditivo <- chmTMB(cbind(success, failure) ~ time + z_ndvi_mean_2120 + (1 | cluster / locationName), data = chm_veado)
ndvi_2560_aditivo <- chmTMB(cbind(success, failure) ~ time + z_ndvi_mean_2560 + (1 | cluster / locationName), data = chm_veado)
ndvi_3000_aditivo <- chmTMB(cbind(success, failure) ~ time + z_ndvi_mean_3000 + (1 | cluster / locationName), data = chm_veado)

ndvi_800_interacao  <- chmTMB(cbind(success, failure) ~ time * z_ndvi_mean_800  + (1 | cluster / locationName), data = chm_veado)
ndvi_1240_interacao <- chmTMB(cbind(success, failure) ~ time * z_ndvi_mean_1240 + (1 | cluster / locationName), data = chm_veado)
ndvi_1680_interacao <- chmTMB(cbind(success, failure) ~ time * z_ndvi_mean_1680 + (1 | cluster / locationName), data = chm_veado)
ndvi_2120_interacao <- chmTMB(cbind(success, failure) ~ time * z_ndvi_mean_2120 + (1 | cluster / locationName), data = chm_veado)
ndvi_2560_interacao <- chmTMB(cbind(success, failure) ~ time * z_ndvi_mean_2560 + (1 | cluster / locationName), data = chm_veado)
ndvi_3000_interacao <- chmTMB(cbind(success, failure) ~ time * z_ndvi_mean_3000 + (1 | cluster / locationName), data = chm_veado)

cat("\n===== B. NDVI =====\n")
cat("Convergencia (TRUE = ok):\n")
print(c(ndvi_800_aditivo  = ndvi_800_aditivo$sdr$pdHess,
        ndvi_1240_aditivo = ndvi_1240_aditivo$sdr$pdHess,
        ndvi_1680_aditivo = ndvi_1680_aditivo$sdr$pdHess,
        ndvi_2120_aditivo = ndvi_2120_aditivo$sdr$pdHess,
        ndvi_2560_aditivo = ndvi_2560_aditivo$sdr$pdHess,
        ndvi_3000_aditivo = ndvi_3000_aditivo$sdr$pdHess,
        ndvi_800_interacao  = ndvi_800_interacao$sdr$pdHess,
        ndvi_1240_interacao = ndvi_1240_interacao$sdr$pdHess,
        ndvi_1680_interacao = ndvi_1680_interacao$sdr$pdHess,
        ndvi_2120_interacao = ndvi_2120_interacao$sdr$pdHess,
        ndvi_2560_interacao = ndvi_2560_interacao$sdr$pdHess,
        ndvi_3000_interacao = ndvi_3000_interacao$sdr$pdHess))

tab_ndvi <- AIC(nulo,
                ndvi_800_aditivo, ndvi_1240_aditivo, ndvi_1680_aditivo,
                ndvi_2120_aditivo, ndvi_2560_aditivo, ndvi_3000_aditivo,
                ndvi_800_interacao, ndvi_1240_interacao, ndvi_1680_interacao,
                ndvi_2120_interacao, ndvi_2560_interacao, ndvi_3000_interacao) %>%
  tibble::rownames_to_column("modelo") %>%
  mutate(variavel = "ndvi_mean",
         n_cameras = n_distinct(chm_veado$locationName))
print(tab_ndvi %>% mutate(dAIC = round(AIC - min(AIC), 2)) %>% arrange(dAIC))


# ============================================================
# C. DISTANCIA ENTRE FRAGMENTOS (enn_mn)
#   enn_mn e NA quando ha menos de 2 fragmentos no buffer. Para o AIC
#   ser comparavel, todas as escalas usam as MESMAS cameras: so as que
#   tem enn_mn em todas as 6 escalas. O nulo e reajustado nesse subconjunto.
# ============================================================
chm_enn <- chm_veado %>%
  filter(!is.na(z_enn_mn_800), !is.na(z_enn_mn_1240), !is.na(z_enn_mn_1680),
         !is.na(z_enn_mn_2120), !is.na(z_enn_mn_2560), !is.na(z_enn_mn_3000))

cat("\n===== C. DISTANCIA ENTRE FRAGMENTOS =====\n")
cat("Cameras usadas:", n_distinct(chm_enn$locationName), "de", n_distinct(chm_veado$locationName),
    "| clusters:", n_distinct(chm_enn$cluster), "\n")
cat("Cameras/sucessos excluidos:\n")
chm_veado %>%
  filter(!locationName %in% chm_enn$locationName) %>%
  summarise(sucessos = sum(success), .by = c(cluster, locationName)) %>%
  print()

nulo_enn <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName), data = chm_enn)

enn_800_aditivo  <- chmTMB(cbind(success, failure) ~ time + z_enn_mn_800  + (1 | cluster / locationName), data = chm_enn)
enn_1240_aditivo <- chmTMB(cbind(success, failure) ~ time + z_enn_mn_1240 + (1 | cluster / locationName), data = chm_enn)
enn_1680_aditivo <- chmTMB(cbind(success, failure) ~ time + z_enn_mn_1680 + (1 | cluster / locationName), data = chm_enn)
enn_2120_aditivo <- chmTMB(cbind(success, failure) ~ time + z_enn_mn_2120 + (1 | cluster / locationName), data = chm_enn)
enn_2560_aditivo <- chmTMB(cbind(success, failure) ~ time + z_enn_mn_2560 + (1 | cluster / locationName), data = chm_enn)
enn_3000_aditivo <- chmTMB(cbind(success, failure) ~ time + z_enn_mn_3000 + (1 | cluster / locationName), data = chm_enn)

enn_800_interacao  <- chmTMB(cbind(success, failure) ~ time * z_enn_mn_800  + (1 | cluster / locationName), data = chm_enn)
enn_1240_interacao <- chmTMB(cbind(success, failure) ~ time * z_enn_mn_1240 + (1 | cluster / locationName), data = chm_enn)
enn_1680_interacao <- chmTMB(cbind(success, failure) ~ time * z_enn_mn_1680 + (1 | cluster / locationName), data = chm_enn)
enn_2120_interacao <- chmTMB(cbind(success, failure) ~ time * z_enn_mn_2120 + (1 | cluster / locationName), data = chm_enn)
enn_2560_interacao <- chmTMB(cbind(success, failure) ~ time * z_enn_mn_2560 + (1 | cluster / locationName), data = chm_enn)
enn_3000_interacao <- chmTMB(cbind(success, failure) ~ time * z_enn_mn_3000 + (1 | cluster / locationName), data = chm_enn)

cat("Convergencia (TRUE = ok):\n")
print(c(nulo_enn = nulo_enn$sdr$pdHess,
        enn_800_aditivo  = enn_800_aditivo$sdr$pdHess,
        enn_1240_aditivo = enn_1240_aditivo$sdr$pdHess,
        enn_1680_aditivo = enn_1680_aditivo$sdr$pdHess,
        enn_2120_aditivo = enn_2120_aditivo$sdr$pdHess,
        enn_2560_aditivo = enn_2560_aditivo$sdr$pdHess,
        enn_3000_aditivo = enn_3000_aditivo$sdr$pdHess,
        enn_800_interacao  = enn_800_interacao$sdr$pdHess,
        enn_1240_interacao = enn_1240_interacao$sdr$pdHess,
        enn_1680_interacao = enn_1680_interacao$sdr$pdHess,
        enn_2120_interacao = enn_2120_interacao$sdr$pdHess,
        enn_2560_interacao = enn_2560_interacao$sdr$pdHess,
        enn_3000_interacao = enn_3000_interacao$sdr$pdHess))

tab_enn <- AIC(nulo_enn,
               enn_800_aditivo, enn_1240_aditivo, enn_1680_aditivo,
               enn_2120_aditivo, enn_2560_aditivo, enn_3000_aditivo,
               enn_800_interacao, enn_1240_interacao, enn_1680_interacao,
               enn_2120_interacao, enn_2560_interacao, enn_3000_interacao) %>%
  tibble::rownames_to_column("modelo") %>%
  mutate(variavel = "enn_mn",
         n_cameras = n_distinct(chm_enn$locationName))
print(tab_enn %>% mutate(dAIC = round(AIC - min(AIC), 2)) %>% arrange(dAIC))


# ============================================================
# D. TABELA FINAL E MELHORES ESCALAS
#   dAIC_nulo = AIC do modelo - AIC do nulo da mesma variavel
#   (negativo = melhor que o nulo; < -2 = suporte claro)
# ============================================================
tab_escala <- bind_rows(tab_native, tab_ndvi, tab_enn) %>%
  mutate(AIC_nulo = AIC[modelo %in% c("nulo", "nulo_enn")][1], .by = variavel) %>%
  mutate(dAIC_nulo = round(AIC - AIC_nulo, 2),
         AIC       = round(AIC, 2),
         escala    = as.numeric(sub("^[a-z]+_([0-9]+)_.*$", "\\1", modelo)),
         efeito    = if_else(modelo %in% c("nulo", "nulo_enn"), "nulo", sub("^.*_", "", modelo))) %>%
  select(variavel, modelo, efeito, escala, df, AIC, dAIC_nulo, n_cameras) %>%
  arrange(variavel, AIC)

write.xlsx(tab_escala, here("output", "tabelas", "03_escala_efeito.xlsx"), rowNames = FALSE)

cat("\n===== D. MELHOR ESCALA DE CADA VARIAVEL =====\n")
melhores <- tab_escala %>%
  filter(!modelo %in% c("nulo", "nulo_enn")) %>%
  slice_min(AIC, n = 1, by = variavel)
print(melhores, row.names = FALSE)
write.xlsx(melhores, here("output", "tabelas", "03_melhores_escalas.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# E. Figura: AIC de cada escala em relacao ao nulo
# ------------------------------------------------------------
figura <- tab_escala %>%
  filter(!modelo %in% c("nulo", "nulo_enn")) %>%
  ggplot(aes(x = escala, y = dAIC_nulo, colour = efeito)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_hline(yintercept = -2, colour = "grey50", linetype = "dashed") +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ variavel, scales = "free_y") +
  scale_x_continuous(breaks = c(800, 1240, 1680, 2120, 2560, 3000)) +
  labs(x = "Raio do buffer (m)",
       y = "AIC (escala) - AIC (nulo)",
       colour = "Efeito",
       caption = "Abaixo da linha tracejada: melhor que o nulo por mais de 2 AIC.") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

ggsave(here("output", "figuras", "03_escala_efeito.png"), figura, width = 9, height = 4, dpi = 300)
print(figura)
