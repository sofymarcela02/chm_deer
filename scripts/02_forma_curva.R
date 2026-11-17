# ============================================================
# 02 - FORMA DA CURVA DE ATIVIDADE (sem covariaveis de paisagem)
# Veado-catingueiro - Iannarilli et al. (2025), Aplicacao 3.1
#
#   A. Estrutura de efeitos aleatorios (cameras aninhadas em clusters)
#   B. Forma da curva: uniforme x unimodal x bimodal
#   C. Campanha (2017 x 2018)
#   D. Figura da curva escolhida
#
# Entrada: data/processed/chm_veado.rds  (gerado por 01_banco_de_dados.R)
# Saidas : output/tabelas/02_*.xlsx, output/figuras/02_forma_curva.png
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

chm_veado <- readRDS(here("data", "processed", "chm_veado.rds"))

# Lembrete do chmTMB:
#   "time" na formula e substituido pelos termos trigonometricos
#   type = "bimodal"  (padrao): cos(t) + sin(t) + cos(2t) + sin(2t)
#   type = "unimodal"         : cos(t) + sin(t)


# ------------------------------------------------------------
# A. Estrutura de efeitos aleatorios (curva bimodal em todos)
#   Cameras aninhadas em clusters: (1 | cluster / locationName)
#   O R expande para (1 | cluster) + (1 | cluster:locationName):
#     - intercepto de cada cluster
#     - intercepto de cada camera dentro do seu cluster
#   No VarCorr os grupos aparecem como "cluster" e "locationName:cluster".
#
#   Intercepto aleatorio  -> frequencia de uso varia entre unidades
#   diag(0 + time | cluster) -> o HORARIO (forma da curva) varia entre
#                               clusters, sem correlacao entre os termos
# ------------------------------------------------------------
re_nenhum <- chmTMB(cbind(success, failure) ~ time,
                    data = chm_veado)

re_camera <- chmTMB(cbind(success, failure) ~ time + (1 | locationName),
                    data = chm_veado)

re_cluster <- chmTMB(cbind(success, failure) ~ time + (1 | cluster),
                     data = chm_veado)

re_aninhado <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                      data = chm_veado)

re_aninhado_horario <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName) +
                                diag(0 + time | cluster),
                              data = chm_veado)

cat("\n===== A. EFEITOS ALEATORIOS =====\n")
cat("Convergencia (TRUE = ok):\n")
print(c(nenhum          = re_nenhum$sdr$pdHess,
        camera          = re_camera$sdr$pdHess,
        cluster         = re_cluster$sdr$pdHess,
        aninhado        = re_aninhado$sdr$pdHess,
        aninhado_horario = re_aninhado_horario$sdr$pdHess))

tab_re <- AIC(re_nenhum, re_camera, re_cluster, re_aninhado, re_aninhado_horario) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_re)
write.xlsx(tibble::rownames_to_column(tab_re, "modelo"), here("output", "tabelas", "02a_efeitos_aleatorios.xlsx"), rowNames = FALSE)

cat("\nDesvios padrao dos efeitos aleatorios (modelo aninhado):\n")
print(VarCorr(re_aninhado))


# ------------------------------------------------------------
# B. Forma da curva
#   Efeitos aleatorios: cameras aninhadas em clusters (desenho amostral)
# ------------------------------------------------------------
uniforme <- chmTMB(cbind(success, failure) ~ 1 + (1 | cluster / locationName),
                   data = chm_veado)

unimodal <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                   type = "unimodal",
                   data = chm_veado)

bimodal <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                  type = "bimodal",
                  data = chm_veado)

cat("\n===== B. FORMA DA CURVA =====\n")
cat("Convergencia:\n")
print(c(uniforme = uniforme$sdr$pdHess,
        unimodal = unimodal$sdr$pdHess,
        bimodal  = bimodal$sdr$pdHess))

tab_forma <- AIC(uniforme, unimodal, bimodal) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_forma)
write.xlsx(tibble::rownames_to_column(tab_forma, "modelo"), here("output", "tabelas", "02b_forma_curva.xlsx"), rowNames = FALSE)

cat("\nTeste de razao de verossimilhanca (uniforme -> unimodal -> bimodal):\n")
print(anova(uniforme, unimodal, bimodal))


# ------------------------------------------------------------
# C. Campanha
#   time + campanha : a campanha muda so o NIVEL de atividade
#   time * campanha : a campanha muda tambem o HORARIO
# ------------------------------------------------------------
bimodal_campanha <- chmTMB(cbind(success, failure) ~ time + campanha +
                             (1 | cluster / locationName),
                           data = chm_veado)

bimodal_x_campanha <- chmTMB(cbind(success, failure) ~ time * campanha +
                               (1 | cluster / locationName),
                             data = chm_veado)

cat("\n===== C. CAMPANHA =====\n")
cat("Convergencia:\n")
print(c(bimodal_campanha   = bimodal_campanha$sdr$pdHess,
        bimodal_x_campanha = bimodal_x_campanha$sdr$pdHess))

tab_campanha <- AIC(bimodal, bimodal_campanha, bimodal_x_campanha) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_campanha)
write.xlsx(tibble::rownames_to_column(tab_campanha, "modelo"), here("output", "tabelas", "02c_campanha.xlsx"), rowNames = FALSE)

cat("\nCoeficientes do modelo bimodal:\n")
print(summary(bimodal)$coefficients$cond)


# ------------------------------------------------------------
# D. Figura: curva bimodal prevista
#
#   predict() do chm_helpers devolve a curva CONDICIONAL: uma camera
#   "tipica", com todos os efeitos aleatorios = 0.
#
#   A curva MARGINAL (media da populacao de cameras) integra sobre os
#   efeitos aleatorios. Aproximacao de Zeger et al. (1988), Biometrics:
#     logit(p_marginal) ~ logit(p_condicional) / sqrt(1 + c^2 * s2)
#     c = 16 * sqrt(3) / (15 * pi);  s2 = soma das variancias aleatorias
#
#   Pontos: proporcao observada em cada hora (todas as cameras juntas)
# ------------------------------------------------------------
var_aleatoria <- VarCorr(bimodal)$cond$cluster[1, 1] +
                 VarCorr(bimodal)$cond$`locationName:cluster`[1, 1]
fator_zeger <- sqrt(1 + (16 * sqrt(3) / (15 * pi))^2 * var_aleatoria)

pred_bimodal <- predict(bimodal, list()) %>%
  mutate(condicional     = response,
         marginal        = plogis(fit / fator_zeger),
         marginal_lcl    = plogis((fit - 1.96 * se.fit) / fator_zeger),
         marginal_ucl    = plogis((fit + 1.96 * se.fit) / fator_zeger))

observado <- chm_veado %>%
  summarise(proporcao = sum(success) / sum(success + failure), .by = timeHour)

cat("\n===== D. ESCALA DAS CURVAS =====\n")
cat("Proporcao observada media      :", signif(mean(observado$proporcao), 3), "\n")
cat("Curva condicional (camera tipica):", signif(mean(pred_bimodal$condicional), 3), "\n")
cat("Curva marginal (Zeger)          :", signif(mean(pred_bimodal$marginal), 3), "\n")

figura <- ggplot(pred_bimodal, aes(x = time)) +
  geom_ribbon(aes(ymin = marginal_lcl, ymax = marginal_ucl), fill = "#1b9e77", alpha = 0.2) +
  geom_line(aes(y = marginal), colour = "#1b9e77", linewidth = 1) +
  geom_point(data = observado, aes(x = timeHour, y = proporcao), colour = "grey30") +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Hora solar",
       y = "Probabilidade de atividade (por hora)",
       title = "Veado-catingueiro - modelo bimodal",
       caption = "Linha: media marginal (aprox. Zeger et al. 1988); faixa: IC 95%. Pontos: proporcao observada.") +
  theme_classic()

figura
ggsave(here("output", "figuras", "02_forma_curva.png"), figura, width = 7, height = 4.5, dpi = 300)
print(figura)
