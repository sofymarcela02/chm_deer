# ============================================================
# 06 - MODELO BINOMIAL INFLADO DE ZEROS (ZIB)
# Veado-catingueiro - 162 cameras - curva bimodal - nativa 1680 m
#
# >>> MODELO FINAL DA ANALISE (decisao de 13/11/2026) <<<
#   zib_interacao: time * z_native_1680 + (1 | cluster / locationName),
#                  ziformula = ~ (1 | locationName)
#   Motivos: usa as 162 cameras; separa "onde" (parte de zeros) de
#   "quando" (parte de atividade); curvas no nivel observado; efeito da
#   nativa robusto sem P150 e P138 (dAIC -4,46).
#   Os scripts 04 (binomial, 162 cameras) e 05 (binomial, so cameras com
#   registro) ficam como analises complementares/de sensibilidade.
#
# Ideia: o modelo tem DUAS partes, ajustadas juntas
#   1. Parte de ZEROS (ziformula): probabilidade de uma camera-hora ser
#      um zero "estrutural" - a camera esta num lugar que o veado nao usa
#      (ou mal posicionada). Com (1 | locationName), essa probabilidade
#      varia entre cameras: separa as cameras "sem veado" das "com veado".
#   2. Parte de ATIVIDADE (formula principal): para as cameras que o veado
#      usa, a probabilidade de atividade em cada hora - a curva bimodal e o
#      efeito da nativa sobre o horario. E a mesma pergunta do script 05,
#      mas sem descartar cameras.
#
# P(sucesso) = (1 - p_zero) x p_atividade
#
#   A. Onde colocar a inflacao de zeros (sem, constante, por camera)
#   B. Com zeros por camera: nulo x aditivo x interacao
#   C. Coeficientes e parte de zeros (quais cameras sao "sem veado")
#   D. Sensibilidade: sem P150, sem P150 e P138
#   E. Comparacao: binomial 162 x binomial so com registro x ZIB
#   F. Figuras
#
# Entrada: data/processed/chm_veado.rds, data/processed/deployments.xlsx
# Saidas : output/tabelas/06_*.xlsx, output/figuras/06_*.png
# Obs.: cada modelo ZIB leva ~10 s; o script todo, alguns minutos.
# ============================================================


# ------------------------------------------------------------
# 0. Pacotes, funcoes do CHM e dados
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(openxlsx)
library(here)
source(here("R", "chm_helpers.R"))

dir.create(here("output", "tabelas"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figuras"), recursive = TRUE, showWarnings = FALSE)

chm_veado   <- readRDS(here("data", "processed", "chm_veado.rds"))
deployments <- read.xlsx(here("data", "processed", "deployments.xlsx"))

# O chmTMB repassa o argumento ziformula para o glmmTMB.
# A palavra "time" so e expandida na formula principal.


# ------------------------------------------------------------
# A. Onde colocar a inflacao de zeros
#   Todos com a parte de atividade = time * z_native_1680
#   zi_nenhum : binomial comum (igual ao script 04)
#   zi_const  : mesma probabilidade de zero estrutural em toda camera-hora
#   zi_nativa : probabilidade de zero estrutural depende da nativa
#   zi_camera : probabilidade de zero estrutural varia entre cameras
#   (zi com nativa + camera foi testado e NAO convergiu)
# ------------------------------------------------------------
zi_nenhum <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                    data = chm_veado)

zi_const <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                   ziformula = ~ 1,
                   data = chm_veado)

zi_nativa <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                    ziformula = ~ z_native_1680,
                    data = chm_veado)

zi_camera <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                    ziformula = ~ (1 | locationName),
                    data = chm_veado)

cat("\n===== A. ESTRUTURA DA INFLACAO DE ZEROS =====\n")
cat("Convergencia (TRUE = ok):\n")
print(c(zi_nenhum = zi_nenhum$sdr$pdHess,
        zi_const  = zi_const$sdr$pdHess,
        zi_nativa = zi_nativa$sdr$pdHess,
        zi_camera = zi_camera$sdr$pdHess))

tab_a <- AIC(zi_nenhum, zi_const, zi_nativa, zi_camera) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_a)
write.xlsx(tibble::rownames_to_column(tab_a, "modelo"),
           here("output", "tabelas", "06a_estrutura_zeros.xlsx"), rowNames = FALSE)

cat("\nDp do efeito de camera na parte de ATIVIDADE:\n")
print(c(zi_nenhum = sqrt(VarCorr(zi_nenhum)$cond$`locationName:cluster`[1, 1]),
        zi_camera = sqrt(VarCorr(zi_camera)$cond$`locationName:cluster`[1, 1])))


# ------------------------------------------------------------
# B. Com zeros por camera: nulo x aditivo x interacao
#   A parte de zeros e a mesma nos 3; muda so a parte de atividade.
# ------------------------------------------------------------
zib_nulo <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                   ziformula = ~ (1 | locationName),
                   data = chm_veado)

zib_aditivo <- chmTMB(cbind(success, failure) ~ time + z_native_1680 + (1 | cluster / locationName),
                      ziformula = ~ (1 | locationName),
                      data = chm_veado)

zib_interacao <- zi_camera

cat("\n===== B. NATIVA (modelo inflado de zeros) =====\n")
cat("Convergencia:\n")
print(c(zib_nulo      = zib_nulo$sdr$pdHess,
        zib_aditivo   = zib_aditivo$sdr$pdHess,
        zib_interacao = zib_interacao$sdr$pdHess))

tab_b <- AIC(zib_nulo, zib_aditivo, zib_interacao) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_b)
write.xlsx(tibble::rownames_to_column(tab_b, "modelo"),
           here("output", "tabelas", "06b_nativa.xlsx"), rowNames = FALSE)

cat("\nLRT (nulo -> aditivo -> interacao):\n")
lrt <- anova(zib_nulo, zib_aditivo, zib_interacao)
print(lrt)
write.xlsx(tibble::rownames_to_column(as.data.frame(lrt), "modelo"),
           here("output", "tabelas", "06b_lrt.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# C. Modelo escolhido: zib_interacao
# ------------------------------------------------------------
cat("\n===== C. MODELO zib_interacao =====\n")
print(summary(zib_interacao))

coef_atividade <- as.data.frame(summary(zib_interacao)$coefficients$cond) %>%
  tibble::rownames_to_column("termo") %>%
  mutate(parte = "atividade") %>%
  relocate(parte)
coef_zeros <- as.data.frame(summary(zib_interacao)$coefficients$zi) %>%
  tibble::rownames_to_column("termo") %>%
  mutate(parte = "zeros") %>%
  relocate(parte)
coeficientes <- bind_rows(coef_atividade, coef_zeros) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
write.xlsx(coeficientes, here("output", "tabelas", "06c_coeficientes.xlsx"), rowNames = FALSE)

# Probabilidade de zero estrutural de cada camera (com o efeito aleatorio
# da camera). Deve ser ~1 nas cameras sem registro e ~0 nas com registro.
# O predict do chm_helpers so faz a curva; para a parte de zeros usamos o
# predict do glmmTMB, por isso a copia com classe "glmmTMB".
zib_glmmTMB <- zib_interacao
class(zib_glmmTMB) <- "glmmTMB"

p_zero_camera <- chm_veado %>%
  mutate(p_zero = predict(zib_glmmTMB, type = "zprob")) %>%
  summarise(p_zero   = mean(p_zero),
            sucessos = sum(success),
            .by = c(cluster, locationName)) %>%
  mutate(com_registro = sucessos > 0)

cat("\nProbabilidade de zero estrutural por camera:\n")
p_zero_camera %>%
  summarise(cameras     = n(),
            p_zero_min  = round(min(p_zero), 3),
            p_zero_med  = round(median(p_zero), 3),
            p_zero_max  = round(max(p_zero), 3),
            .by = com_registro) %>%
  print()
write.xlsx(p_zero_camera, here("output", "tabelas", "06c_prob_zero_por_camera.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# D. Sensibilidade: P150 (97,5% nativa) e P138 (67,2% nativa)
# ------------------------------------------------------------
cat("\n===== D. SENSIBILIDADE =====\n")

chm_sem_P150 <- chm_veado %>% filter(cluster != "P150")
zib_nulo_sem_P150 <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                            ziformula = ~ (1 | locationName),
                            data = chm_sem_P150)
zib_interacao_sem_P150 <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                                 ziformula = ~ (1 | locationName),
                                 data = chm_sem_P150)

chm_sem_P150_P138 <- chm_veado %>% filter(!cluster %in% c("P150", "P138"))
zib_nulo_sem_P150_P138 <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                                 ziformula = ~ (1 | locationName),
                                 data = chm_sem_P150_P138)
zib_interacao_sem_P150_P138 <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                                      ziformula = ~ (1 | locationName),
                                      data = chm_sem_P150_P138)

cat("Convergencia:\n")
print(c(zib_nulo_sem_P150           = zib_nulo_sem_P150$sdr$pdHess,
        zib_interacao_sem_P150      = zib_interacao_sem_P150$sdr$pdHess,
        zib_nulo_sem_P150_P138      = zib_nulo_sem_P150_P138$sdr$pdHess,
        zib_interacao_sem_P150_P138 = zib_interacao_sem_P150_P138$sdr$pdHess))

tab_d <- data.frame(
  cenario  = c("todos os clusters", "sem P150", "sem P150 e P138"),
  clusters = c(n_distinct(chm_veado$cluster),
               n_distinct(chm_sem_P150$cluster),
               n_distinct(chm_sem_P150_P138$cluster)),
  sucessos = c(sum(chm_veado$success),
               sum(chm_sem_P150$success),
               sum(chm_sem_P150_P138$success)),
  AIC_nulo      = c(AIC(zib_nulo), AIC(zib_nulo_sem_P150), AIC(zib_nulo_sem_P150_P138)),
  AIC_interacao = c(AIC(zib_interacao), AIC(zib_interacao_sem_P150), AIC(zib_interacao_sem_P150_P138))
) %>%
  mutate(dAIC_interacao_nulo = round(AIC_interacao - AIC_nulo, 2),
         AIC_nulo            = round(AIC_nulo, 2),
         AIC_interacao       = round(AIC_interacao, 2))
print(tab_d, row.names = FALSE)
write.xlsx(tab_d, here("output", "tabelas", "06d_sensibilidade.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# E. Comparacao dos 3 caminhos
#   binomial com 162 cameras (script 04)
#   binomial so com as cameras com registro (script 05)
#   binomial inflado de zeros com 162 cameras (este script)
# ------------------------------------------------------------
cat("\n===== E. COMPARACAO =====\n")

cameras_com_registro <- chm_veado %>%
  summarise(sucessos = sum(success), .by = locationName) %>%
  filter(sucessos > 0)
chm_reg <- chm_veado %>% filter(locationName %in% cameras_com_registro$locationName)

bin_reg_nulo <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                       data = chm_reg)
bin_reg_interacao <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                            data = chm_reg)
bin_162_nulo <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                       data = chm_veado)

tab_e <- data.frame(
  abordagem = c("binomial - 162 cameras", "binomial - so com registro", "ZIB - 162 cameras"),
  cameras   = c(162, n_distinct(chm_reg$locationName), 162),
  dAIC_interacao_nulo = c(AIC(zi_nenhum) - AIC(bin_162_nulo),
                          AIC(bin_reg_interacao) - AIC(bin_reg_nulo),
                          AIC(zib_interacao) - AIC(zib_nulo)),
  dp_camera_atividade = c(sqrt(VarCorr(zi_nenhum)$cond$`locationName:cluster`[1, 1]),
                          sqrt(VarCorr(bin_reg_interacao)$cond$`locationName:cluster`[1, 1]),
                          sqrt(VarCorr(zib_interacao)$cond$`locationName:cluster`[1, 1]))
) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
print(tab_e, row.names = FALSE)
write.xlsx(tab_e, here("output", "tabelas", "06e_comparacao_abordagens.xlsx"), rowNames = FALSE)

cat("\nCoeficientes da parte de atividade lado a lado:\n")
tab_coef <- data.frame(
  termo           = names(fixef(zib_interacao)$cond),
  binomial_162    = round(fixef(zi_nenhum)$cond, 3),
  binomial_reg    = round(fixef(bin_reg_interacao)$cond, 3),
  ZIB_162         = round(fixef(zib_interacao)$cond, 3),
  row.names       = NULL
)
print(tab_coef, row.names = FALSE)
write.xlsx(tab_coef, here("output", "tabelas", "06e_coeficientes_lado_a_lado.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# F. Figuras - parte de ATIVIDADE do ZIB
#   Curva de atividade para cameras que o veado USA (parte de zeros
#   fora). Marginal sobre os efeitos aleatorios da parte de atividade
#   (aprox. de Zeger et al. 1988; aqui o efeito aleatorio e pequeno).
#   Pontos: proporcao observada so nas cameras com registro.
# ------------------------------------------------------------
var_aleatoria <- VarCorr(zib_interacao)$cond$cluster[1, 1] +
                 VarCorr(zib_interacao)$cond$`locationName:cluster`[1, 1]
fator_zeger <- sqrt(1 + (16 * sqrt(3) / (15 * pi))^2 * var_aleatoria)

nativa_55  <- deployments %>% distinct(cluster, native_1680, z_native_1680)
niveis_z   <- quantile(nativa_55$z_native_1680, c(0.10, 0.50, 0.90))
niveis_pct <- round(quantile(nativa_55$native_1680, c(0.10, 0.50, 0.90)))

rotulo_pouca <- paste0("Pouca (", niveis_pct[1], "% nativa)")
rotulo_media <- paste0("Media (", niveis_pct[2], "% nativa)")
rotulo_muita <- paste0("Muita (", niveis_pct[3], "% nativa)")

# predict do chm_helpers devolve a parte de atividade (escala logit)
pred <- predict(zib_interacao, list(z_native_1680 = niveis_z)) %>%
  mutate(nivel = case_when(z_native_1680 == niveis_z[1] ~ rotulo_pouca,
                           z_native_1680 == niveis_z[2] ~ rotulo_media,
                           z_native_1680 == niveis_z[3] ~ rotulo_muita),
         nivel = factor(nivel, levels = c(rotulo_pouca, rotulo_media, rotulo_muita)),
         marginal     = plogis(fit / fator_zeger),
         marginal_lcl = plogis((fit - 1.96 * se.fit) / fator_zeger),
         marginal_ucl = plogis((fit + 1.96 * se.fit) / fator_zeger)) %>%
  mutate(relativo = marginal / sum(marginal), .by = nivel)

cores <- c("#d95f02", "#7f7f7f", "#1b9e77")

figura_prob <- ggplot(pred, aes(x = time, colour = nivel, fill = nivel)) +
  geom_ribbon(aes(ymin = marginal_lcl, ymax = marginal_ucl), colour = NA, alpha = 0.15) +
  geom_line(aes(y = marginal), linewidth = 1) +
  facet_wrap(~ nivel, nrow = 1) +
  scale_colour_manual(values = cores, aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Hora solar", y = "Probabilidade de atividade (por hora)",
       title = "Veado-catingueiro - binomial inflado de zeros - vegetacao nativa em 1680 m",
       caption = "Parte de atividade (cameras usadas pelo veado); media marginal; faixa: IC 95%.") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(here("output", "figuras", "06_zib_atividade_nativa_probabilidade.png"), figura_prob,
       width = 10, height = 4, dpi = 300)

figura_rel <- ggplot(pred, aes(x = time, y = relativo, colour = nivel)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = cores) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Hora solar", y = "Proporcao da atividade diaria", colour = NULL,
       title = "Distribuicao da atividade ao longo do dia - binomial inflado de zeros",
       caption = "Cada curva soma 1 ao longo das 24 h (so a forma).") +
  theme_classic() +
  theme(legend.position = "top")
ggsave(here("output", "figuras", "06_zib_atividade_nativa_relativa.png"), figura_rel,
       width = 7, height = 4.5, dpi = 300)

print(figura_prob)
print(figura_rel)

cat("\nNivel medio da curva (nivel medio de nativa) x observado nas cameras com registro:\n")
cat("Curva ZIB (Media):", signif(mean(pred$marginal[pred$nivel == rotulo_media]), 3), "\n")
cat("Observado        :", signif(sum(chm_reg$success) / sum(chm_reg$success + chm_reg$failure), 3), "\n")

cat("\nHorario do pico principal (hora solar):\n")
pred %>%
  slice_max(marginal, n = 1, by = nivel) %>%
  transmute(nivel, hora_pico = round(time, 1)) %>%
  print(row.names = FALSE)
