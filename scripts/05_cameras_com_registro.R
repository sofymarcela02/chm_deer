# ============================================================
# 05 - ANALISE SO COM AS CAMERAS QUE REGISTRARAM O VEADO
#
# >>> ANALISE COMPLEMENTAR. O modelo principal e o binomial inflado de
#     zeros do script 06 (decisao de 13/11/2026), que usa as 162 cameras.
# Veado-catingueiro - curva bimodal, cameras aninhadas em clusters
#
# Por que: 134 das 162 cameras nunca registraram o veado. Essas cameras
# nao informam sobre o HORARIO da atividade, mas inflam o efeito
# aleatorio de camera (dp ~ 4,6) e distorcem o nivel das curvas.
# Aqui a pergunta passa a ser: "QUANDO o veado e ativo, ONDE ele ocorre?"
#
#   A. Estrutura de efeitos aleatorios
#   B. Forma da curva (uniforme x unimodal x bimodal)
#   C. Vegetacao nativa 1680 m: nulo x aditivo x interacao
#   D. R2 marginal e condicional (Nakagawa)
#   E. Sensibilidade: sem P150, sem P150 e P138
#   F. Comparacao com a analise com todas as 162 cameras (script 04)
#   G. Figuras
#
# Obs.: as covariaveis z_ sao as mesmas do script 01 (padronizadas com
# os 55 clusters), para os coeficientes serem comparaveis com o script 04.
#
# Entrada: data/processed/chm_veado.rds, data/processed/deployments.xlsx
# Saidas : output/tabelas/05_*.xlsx, output/figuras/05_*.png
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

# Cameras com pelo menos 1 sucesso (1 dia com registro de veado)
cameras_com_registro <- chm_veado %>%
  summarise(sucessos = sum(success), .by = c(cluster, locationName)) %>%
  filter(sucessos > 0)

chm_reg <- chm_veado %>%
  filter(locationName %in% cameras_com_registro$locationName)

cat("\n===== DADOS =====\n")
cat("Cameras com registro:", n_distinct(chm_reg$locationName), "de", n_distinct(chm_veado$locationName), "\n")
cat("Clusters            :", n_distinct(chm_reg$cluster), "de", n_distinct(chm_veado$cluster), "\n")
cat("Linhas da base      :", nrow(chm_reg), "\n")
cat("Sucessos            :", sum(chm_reg$success), "(todos os", sum(chm_veado$success), "sucessos ficam)\n")
cat("Ocasioes            :", sum(chm_reg$success + chm_reg$failure), "\n")

cat("\nSucessos por camera:\n")
print(table(sucessos_por_camera = cameras_com_registro$sucessos))

# Vegetacao nativa dos clusters com registro
nativa_reg <- chm_reg %>% distinct(cluster, native_1680, z_native_1680)
cat("\nVegetacao nativa 1680 m nos", nrow(nativa_reg), "clusters com registro:\n")
print(summary(nativa_reg$native_1680))


# ------------------------------------------------------------
# A. Estrutura de efeitos aleatorios (curva bimodal)
# ------------------------------------------------------------
re_nenhum <- chmTMB(cbind(success, failure) ~ time,
                    data = chm_reg)

re_camera <- chmTMB(cbind(success, failure) ~ time + (1 | locationName),
                    data = chm_reg)

re_cluster <- chmTMB(cbind(success, failure) ~ time + (1 | cluster),
                     data = chm_reg)

re_aninhado <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                      data = chm_reg)

cat("\n===== A. EFEITOS ALEATORIOS =====\n")
cat("Convergencia (TRUE = ok):\n")
print(c(nenhum   = re_nenhum$sdr$pdHess,
        camera   = re_camera$sdr$pdHess,
        cluster  = re_cluster$sdr$pdHess,
        aninhado = re_aninhado$sdr$pdHess))

tab_a <- AIC(re_nenhum, re_camera, re_cluster, re_aninhado) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_a)
write.xlsx(tibble::rownames_to_column(tab_a, "modelo"),
           here("output", "tabelas", "05a_efeitos_aleatorios.xlsx"), rowNames = FALSE)

cat("\nDesvios padrao dos efeitos aleatorios (modelo aninhado):\n")
print(VarCorr(re_aninhado))


# ------------------------------------------------------------
# B. Forma da curva (cameras aninhadas em clusters)
# ------------------------------------------------------------
uniforme <- chmTMB(cbind(success, failure) ~ 1 + (1 | cluster / locationName),
                   data = chm_reg)

unimodal <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                   type = "unimodal",
                   data = chm_reg)

bimodal <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                  type = "bimodal",
                  data = chm_reg)

cat("\n===== B. FORMA DA CURVA =====\n")
cat("Convergencia:\n")
print(c(uniforme = uniforme$sdr$pdHess,
        unimodal = unimodal$sdr$pdHess,
        bimodal  = bimodal$sdr$pdHess))

tab_b <- AIC(uniforme, unimodal, bimodal) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_b)
write.xlsx(tibble::rownames_to_column(tab_b, "modelo"),
           here("output", "tabelas", "05b_forma_curva.xlsx"), rowNames = FALSE)

cat("\nLRT (uniforme -> unimodal -> bimodal):\n")
print(anova(uniforme, unimodal, bimodal))


# ------------------------------------------------------------
# C. Vegetacao nativa (1680 m)
#   nulo      = bimodal (bloco B)
#   aditivo   : a nativa muda so o NIVEL de atividade
#   interacao : a nativa muda tambem o HORARIO
# ------------------------------------------------------------
nulo <- bimodal

nativa_aditivo <- chmTMB(cbind(success, failure) ~ time + z_native_1680 + (1 | cluster / locationName),
                         data = chm_reg)

nativa_interacao <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                           data = chm_reg)

cat("\n===== C. VEGETACAO NATIVA =====\n")
cat("Convergencia:\n")
print(c(nativa_aditivo   = nativa_aditivo$sdr$pdHess,
        nativa_interacao = nativa_interacao$sdr$pdHess))

tab_c <- AIC(nulo, nativa_aditivo, nativa_interacao) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_c)
write.xlsx(tibble::rownames_to_column(tab_c, "modelo"),
           here("output", "tabelas", "05c_nativa.xlsx"), rowNames = FALSE)

cat("\nLRT (nulo -> aditivo -> interacao):\n")
lrt <- anova(nulo, nativa_aditivo, nativa_interacao)
print(lrt)
write.xlsx(tibble::rownames_to_column(as.data.frame(lrt), "modelo"),
           here("output", "tabelas", "05c_lrt.xlsx"), rowNames = FALSE)

cat("\nModelo nativa_interacao:\n")
print(summary(nativa_interacao))

coeficientes <- as.data.frame(summary(nativa_interacao)$coefficients$cond) %>%
  tibble::rownames_to_column("termo") %>%
  mutate(across(-termo, ~ signif(.x, 3)))
write.xlsx(coeficientes, here("output", "tabelas", "05c_coeficientes.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# D. R2 marginal e condicional (Nakagawa & Schielzeth 2013)
#   Mesmo calculo do script 04 (escala logit, residual = pi^2/3)
# ------------------------------------------------------------
cat("\n===== D. R2 DE NAKAGAWA =====\n")

var_fixos_nulo     <- var(as.vector(getME(nulo, "X") %*% fixef(nulo)$cond))
var_aleatoria_nulo <- VarCorr(nulo)$cond$cluster[1, 1] +
                      VarCorr(nulo)$cond$`locationName:cluster`[1, 1]

var_fixos_aditivo     <- var(as.vector(getME(nativa_aditivo, "X") %*% fixef(nativa_aditivo)$cond))
var_aleatoria_aditivo <- VarCorr(nativa_aditivo)$cond$cluster[1, 1] +
                         VarCorr(nativa_aditivo)$cond$`locationName:cluster`[1, 1]

var_fixos_interacao     <- var(as.vector(getME(nativa_interacao, "X") %*% fixef(nativa_interacao)$cond))
var_aleatoria_interacao <- VarCorr(nativa_interacao)$cond$cluster[1, 1] +
                           VarCorr(nativa_interacao)$cond$`locationName:cluster`[1, 1]

var_residual <- pi^2 / 3

tab_d <- data.frame(
  modelo        = c("nulo", "nativa_aditivo", "nativa_interacao"),
  var_fixos     = c(var_fixos_nulo, var_fixos_aditivo, var_fixos_interacao),
  var_aleatoria = c(var_aleatoria_nulo, var_aleatoria_aditivo, var_aleatoria_interacao),
  var_residual  = var_residual
) %>%
  mutate(R2_marginal    = var_fixos / (var_fixos + var_aleatoria + var_residual),
         R2_condicional = (var_fixos + var_aleatoria) / (var_fixos + var_aleatoria + var_residual),
         across(-modelo, ~ round(.x, 3)))
print(tab_d, row.names = FALSE)
write.xlsx(tab_d, here("output", "tabelas", "05d_r2_nakagawa.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# E. Sensibilidade: P150 (97,5% nativa) e P138 (67,2% nativa)
#   Os 2 pontos com muita nativa e com registro (ver script 04)
# ------------------------------------------------------------
cat("\n===== E. SENSIBILIDADE =====\n")

chm_reg_sem_P150 <- chm_reg %>% filter(cluster != "P150")
nulo_sem_P150 <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                        data = chm_reg_sem_P150)
interacao_sem_P150 <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                             data = chm_reg_sem_P150)

chm_reg_sem_P150_P138 <- chm_reg %>% filter(!cluster %in% c("P150", "P138"))
nulo_sem_P150_P138 <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                             data = chm_reg_sem_P150_P138)
interacao_sem_P150_P138 <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                                  data = chm_reg_sem_P150_P138)

cat("Convergencia:\n")
print(c(nulo_sem_P150           = nulo_sem_P150$sdr$pdHess,
        interacao_sem_P150      = interacao_sem_P150$sdr$pdHess,
        nulo_sem_P150_P138      = nulo_sem_P150_P138$sdr$pdHess,
        interacao_sem_P150_P138 = interacao_sem_P150_P138$sdr$pdHess))

tab_e <- data.frame(
  cenario  = c("cameras com registro", "sem P150", "sem P150 e P138"),
  cameras  = c(n_distinct(chm_reg$locationName),
               n_distinct(chm_reg_sem_P150$locationName),
               n_distinct(chm_reg_sem_P150_P138$locationName)),
  clusters = c(n_distinct(chm_reg$cluster),
               n_distinct(chm_reg_sem_P150$cluster),
               n_distinct(chm_reg_sem_P150_P138$cluster)),
  sucessos = c(sum(chm_reg$success),
               sum(chm_reg_sem_P150$success),
               sum(chm_reg_sem_P150_P138$success)),
  AIC_nulo      = c(AIC(nulo), AIC(nulo_sem_P150), AIC(nulo_sem_P150_P138)),
  AIC_interacao = c(AIC(nativa_interacao), AIC(interacao_sem_P150), AIC(interacao_sem_P150_P138))
) %>%
  mutate(dAIC_interacao_nulo = round(AIC_interacao - AIC_nulo, 2),
         AIC_nulo            = round(AIC_nulo, 2),
         AIC_interacao       = round(AIC_interacao, 2))
print(tab_e, row.names = FALSE)
write.xlsx(tab_e, here("output", "tabelas", "05e_sensibilidade.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# F. Comparacao: 162 cameras (script 04) x cameras com registro
#   O AIC nao se compara entre bases diferentes; compara-se o dAIC
#   (interacao - nulo) dentro de cada base, o dp dos efeitos aleatorios,
#   o R2 e o nivel das curvas.
# ------------------------------------------------------------
cat("\n===== F. COMPARACAO COM AS 162 CAMERAS =====\n")

nulo_162 <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                   data = chm_veado)
interacao_162 <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                        data = chm_veado)

var_fixos_162     <- var(as.vector(getME(interacao_162, "X") %*% fixef(interacao_162)$cond))
var_aleatoria_162 <- VarCorr(interacao_162)$cond$cluster[1, 1] +
                     VarCorr(interacao_162)$cond$`locationName:cluster`[1, 1]

tab_f <- data.frame(
  base      = c("162 cameras (todas)", "cameras com registro"),
  cameras   = c(n_distinct(chm_veado$locationName), n_distinct(chm_reg$locationName)),
  clusters  = c(n_distinct(chm_veado$cluster), n_distinct(chm_reg$cluster)),
  dAIC_interacao_nulo = c(AIC(interacao_162) - AIC(nulo_162),
                          AIC(nativa_interacao) - AIC(nulo)),
  dp_camera  = c(sqrt(VarCorr(interacao_162)$cond$`locationName:cluster`[1, 1]),
                 sqrt(VarCorr(nativa_interacao)$cond$`locationName:cluster`[1, 1])),
  dp_cluster = c(sqrt(VarCorr(interacao_162)$cond$cluster[1, 1]),
                 sqrt(VarCorr(nativa_interacao)$cond$cluster[1, 1])),
  R2_marginal = c(var_fixos_162 / (var_fixos_162 + var_aleatoria_162 + var_residual),
                  var_fixos_interacao / (var_fixos_interacao + var_aleatoria_interacao + var_residual)),
  proporcao_observada = c(sum(chm_veado$success) / sum(chm_veado$success + chm_veado$failure),
                          sum(chm_reg$success) / sum(chm_reg$success + chm_reg$failure))
) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
print(tab_f, row.names = FALSE)
write.xlsx(tab_f, here("output", "tabelas", "05f_comparacao_162_x_registro.xlsx"), rowNames = FALSE)

cat("\nCoeficientes de interacao lado a lado:\n")
data.frame(
  termo        = names(fixef(nativa_interacao)$cond),
  todas_162    = round(fixef(interacao_162)$cond, 3),
  com_registro = round(fixef(nativa_interacao)$cond, 3),
  row.names    = NULL
) %>% print(row.names = FALSE)


# ------------------------------------------------------------
# G. Figuras
#   Curva MARGINAL (aprox. de Zeger et al. 1988). Com as cameras com
#   registro o efeito aleatorio e pequeno, entao a curva marginal fica
#   proxima da condicional e no nivel observado.
# ------------------------------------------------------------

# G1. Curva geral (modelo bimodal sem nativa) x proporcao observada
var_aleatoria_bimodal <- VarCorr(bimodal)$cond$cluster[1, 1] +
                         VarCorr(bimodal)$cond$`locationName:cluster`[1, 1]
fator_zeger_bimodal <- sqrt(1 + (16 * sqrt(3) / (15 * pi))^2 * var_aleatoria_bimodal)

pred_geral <- predict(bimodal, list()) %>%
  mutate(marginal     = plogis(fit / fator_zeger_bimodal),
         marginal_lcl = plogis((fit - 1.96 * se.fit) / fator_zeger_bimodal),
         marginal_ucl = plogis((fit + 1.96 * se.fit) / fator_zeger_bimodal))

observado <- chm_reg %>%
  summarise(proporcao = sum(success) / sum(success + failure), .by = timeHour)

cat("\n===== G. NIVEL DAS CURVAS =====\n")
cat("Proporcao observada media:", signif(mean(observado$proporcao), 3), "\n")
cat("Curva marginal media     :", signif(mean(pred_geral$marginal), 3), "\n")

figura_geral <- ggplot(pred_geral, aes(x = time)) +
  geom_ribbon(aes(ymin = marginal_lcl, ymax = marginal_ucl), fill = "#1b9e77", alpha = 0.2) +
  geom_line(aes(y = marginal), colour = "#1b9e77", linewidth = 1) +
  geom_point(data = observado, aes(x = timeHour, y = proporcao), colour = "grey30") +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Hora solar", y = "Probabilidade de atividade (por hora)",
       title = "Veado-catingueiro - cameras com registro (modelo bimodal)",
       caption = "Linha: media marginal; faixa: IC 95%. Pontos: proporcao observada.") +
  theme_classic()
ggsave(here("output", "figuras", "05_curva_geral.png"), figura_geral, width = 7, height = 4.5, dpi = 300)
print(figura_geral)

# G2. Atividade com pouca, media e muita vegetacao nativa
#   Mesmos niveis do script 04 (percentis 10%, 50%, 90% dos 55 clusters)
nativa_55 <- deployments %>% distinct(cluster, native_1680, z_native_1680)
niveis_z   <- quantile(nativa_55$z_native_1680, c(0.10, 0.50, 0.90))
niveis_pct <- round(quantile(nativa_55$native_1680, c(0.10, 0.50, 0.90)))

fator_zeger <- sqrt(1 + (16 * sqrt(3) / (15 * pi))^2 * var_aleatoria_interacao)

rotulo_pouca <- paste0("Pouca (", niveis_pct[1], "% nativa)")
rotulo_media <- paste0("Media (", niveis_pct[2], "% nativa)")
rotulo_muita <- paste0("Muita (", niveis_pct[3], "% nativa)")

pred_nativa <- predict(nativa_interacao, list(z_native_1680 = niveis_z)) %>%
  mutate(nivel = case_when(z_native_1680 == niveis_z[1] ~ rotulo_pouca,
                           z_native_1680 == niveis_z[2] ~ rotulo_media,
                           z_native_1680 == niveis_z[3] ~ rotulo_muita),
         nivel = factor(nivel, levels = c(rotulo_pouca, rotulo_media, rotulo_muita)),
         marginal     = plogis(fit / fator_zeger),
         marginal_lcl = plogis((fit - 1.96 * se.fit) / fator_zeger),
         marginal_ucl = plogis((fit + 1.96 * se.fit) / fator_zeger)) %>%
  mutate(relativo = marginal / sum(marginal), .by = nivel)

cores <- c("#d95f02", "#7f7f7f", "#1b9e77")

figura_prob <- ggplot(pred_nativa, aes(x = time, colour = nivel, fill = nivel)) +
  geom_ribbon(aes(ymin = marginal_lcl, ymax = marginal_ucl), colour = NA, alpha = 0.15) +
  geom_line(aes(y = marginal), linewidth = 1) +
  facet_wrap(~ nivel, nrow = 1) +
  scale_colour_manual(values = cores, aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Hora solar", y = "Probabilidade de atividade (por hora)",
       title = "Veado-catingueiro - cameras com registro - vegetacao nativa em 1680 m",
       caption = "Media marginal (aprox. Zeger et al. 1988); faixa: IC 95%.") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(here("output", "figuras", "05_atividade_nativa_probabilidade.png"), figura_prob,
       width = 10, height = 4, dpi = 300)

figura_rel <- ggplot(pred_nativa, aes(x = time, y = relativo, colour = nivel)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = cores) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Hora solar", y = "Proporcao da atividade diaria", colour = NULL,
       title = "Distribuicao da atividade ao longo do dia - cameras com registro",
       caption = "Cada curva soma 1 ao longo das 24 h (so a forma).") +
  theme_classic() +
  theme(legend.position = "top")
ggsave(here("output", "figuras", "05_atividade_nativa_relativa.png"), figura_rel,
       width = 7, height = 4.5, dpi = 300)

print(figura_prob)
print(figura_rel)

cat("\nHorario do pico principal (hora solar):\n")
pred_nativa %>%
  slice_max(marginal, n = 1, by = nivel) %>%
  transmute(nivel, hora_pico = round(time, 1)) %>%
  print(row.names = FALSE)
