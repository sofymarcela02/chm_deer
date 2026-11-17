# ============================================================
# 04 - MODELO COM VEGETACAO NATIVA (1680 m) - BINOMIAL, 162 CAMERAS
#
# >>> ANALISE COMPLEMENTAR. O modelo principal e o binomial inflado de
#     zeros do script 06 (decisao de 13/11/2026). Aqui o excesso de
#     cameras sem registro infla o efeito aleatorio de camera (dp ~4,6).
# Veado-catingueiro - curva bimodal, cameras aninhadas em clusters
#
#   A. Nulo x aditivo x interacao (AIC e LRT)
#   B. Horario variando entre clusters (inclinacao aleatoria)
#   C. Coeficientes e efeitos aleatorios do modelo escolhido
#   C2. R2 marginal e condicional (Nakagawa)
#   D. Sensibilidade: sem P150, sem P150 e P138, nativa transformada
#   E. Figuras: atividade com pouca, media e muita vegetacao nativa
#
# Entrada: data/processed/chm_veado.rds, data/processed/deployments.xlsx
# Saidas : output/tabelas/04_*.xlsx, output/figuras/04_*.png
#
# ATENCAO - PONTOS P150 E P138
#   O efeito da nativa depende destes 2 pontos (ver bloco D).
#   - P150: 97,5% de nativa (z = 3,54), o MAIOR valor do estudo e um
#           extremo (mediana dos pontos = 17%). Campanha 2018; cerrado/cerradao.
#           5 sucessos (C2: 1; C3: 4), as 8:00-10:00 e 18:00-19:00.
#   - P138: 67,2% de nativa (z = 2,05), 4o maior valor. Campanha 2018;
#           transicao. 4 sucessos, todos na C1, entre 7:00 e 13:00.
#   Sao os unicos pontos com MUITA nativa que tem registro de veado:
#   P148 (74%) e P147 (68%) tem mais nativa que o P138, mas nenhum registro.
#   Todos os 9 sucessos sao diurnos/crepusculares; e isso que gera os
#   "dois picos" previstos para paisagens com muita nativa.
#   Os dois tambem nao tem enn_mn (menos de 2 fragmentos no buffer), por
#   isso saem de qualquer analise que inclua a distancia entre fragmentos.
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

# Valores de vegetacao nativa por cluster (1 linha por ponto)
nativa_cluster <- deployments %>%
  distinct(cluster, native_1680, z_native_1680)

cat("\nVegetacao nativa em 1680 m (55 clusters):\n")
print(summary(nativa_cluster$native_1680))


# ------------------------------------------------------------
# A. Nulo x aditivo x interacao
#   nulo      : so a curva bimodal
#   aditivo   : a nativa muda so o NIVEL de atividade
#   interacao : a nativa muda tambem o HORARIO (forma da curva)
# ------------------------------------------------------------
nulo <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
               data = chm_veado)

nativa_aditivo <- chmTMB(cbind(success, failure) ~ time + z_native_1680 + (1 | cluster / locationName),
                         data = chm_veado)

nativa_interacao <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                           data = chm_veado)

cat("\n===== A. NULO x ADITIVO x INTERACAO =====\n")
cat("Convergencia (TRUE = ok):\n")
print(c(nulo             = nulo$sdr$pdHess,
        nativa_aditivo   = nativa_aditivo$sdr$pdHess,
        nativa_interacao = nativa_interacao$sdr$pdHess))

tab_a <- AIC(nulo, nativa_aditivo, nativa_interacao) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_a)
write.xlsx(tibble::rownames_to_column(tab_a, "modelo"), here("output", "tabelas", "04a_nulo_aditivo_interacao.xlsx"), rowNames = FALSE)

cat("\nTeste de razao de verossimilhanca (nulo -> aditivo -> interacao):\n")
lrt <- anova(nulo, nativa_aditivo, nativa_interacao)
print(lrt)
write.xlsx(tibble::rownames_to_column(as.data.frame(lrt), "modelo"), here("output", "tabelas", "04a_lrt.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# B. Horario variando entre clusters
#   diag(0 + time | cluster): cada cluster pode ter a curva deslocada,
#   alem do efeito da nativa (termos sem correlacao entre si)
# ------------------------------------------------------------
nativa_interacao_horario <- chmTMB(cbind(success, failure) ~ time * z_native_1680 +
                                     (1 | cluster / locationName) + diag(0 + time | cluster),
                                   data = chm_veado)

cat("\n===== B. INCLINACAO ALEATORIA NO HORARIO =====\n")
cat("Convergencia:", nativa_interacao_horario$sdr$pdHess, "\n")
tab_b <- AIC(nativa_interacao, nativa_interacao_horario) %>%
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC  = round(AIC, 2)) %>%
  arrange(dAIC)
print(tab_b)
write.xlsx(tibble::rownames_to_column(tab_b, "modelo"), here("output", "tabelas", "04b_inclinacao_aleatoria.xlsx"), rowNames = FALSE)


# ------------------------------------------------------------
# C. Modelo escolhido: nativa_interacao
# ------------------------------------------------------------
cat("\n===== C. MODELO ESCOLHIDO: nativa_interacao =====\n")
print(summary(nativa_interacao))

coeficientes <- as.data.frame(summary(nativa_interacao)$coefficients$cond) %>%
  tibble::rownames_to_column("termo") %>%
  mutate(across(-termo, ~ signif(.x, 3)))
write.xlsx(coeficientes, here("output", "tabelas", "04c_coeficientes.xlsx"), rowNames = FALSE)

# Teste conjunto dos 4 termos de interacao (Wald), complementa o LRT
termos_interacao <- grep(":", names(fixef(nativa_interacao)$cond), value = TRUE)
b_int <- fixef(nativa_interacao)$cond[termos_interacao]
V_int <- vcov(nativa_interacao)$cond[termos_interacao, termos_interacao]
wald  <- as.numeric(t(b_int) %*% solve(V_int) %*% b_int)
cat("\nWald conjunto das interacoes: X2 =", round(wald, 2),
    "| gl =", length(b_int),
    "| p =", signif(pchisq(wald, length(b_int), lower.tail = FALSE), 3), "\n")


# ------------------------------------------------------------
# C2. R2 marginal e condicional (Nakagawa & Schielzeth 2013)
#   Escala latente (logit). Componentes de variancia:
#     var_fixos     = variancia do preditor linear dos efeitos fixos
#                     (X %*% beta) entre as linhas da base
#     var_aleatoria = cluster + camera dentro do cluster
#     var_residual  = pi^2 / 3  (variancia da distribuicao logistica)
#   R2 marginal    = var_fixos / (var_fixos + var_aleatoria + pi^2/3)
#                    -> o quanto os EFEITOS FIXOS explicam
#   R2 condicional = (var_fixos + var_aleatoria) / (mesmo denominador)
#                    -> o quanto fixos + aleatorios explicam juntos
#   No nulo, os efeitos fixos sao so a curva de atividade (time); a
#   diferenca entre os modelos mostra o que a nativa acrescenta.
# ------------------------------------------------------------
cat("\n===== C2. R2 DE NAKAGAWA =====\n")

# Nulo
var_fixos_nulo     <- var(as.vector(getME(nulo, "X") %*% fixef(nulo)$cond))
var_aleatoria_nulo <- VarCorr(nulo)$cond$cluster[1, 1] +
                      VarCorr(nulo)$cond$`locationName:cluster`[1, 1]

# Aditivo
var_fixos_aditivo     <- var(as.vector(getME(nativa_aditivo, "X") %*% fixef(nativa_aditivo)$cond))
var_aleatoria_aditivo <- VarCorr(nativa_aditivo)$cond$cluster[1, 1] +
                         VarCorr(nativa_aditivo)$cond$`locationName:cluster`[1, 1]

# Interacao (modelo escolhido)
var_fixos_interacao     <- var(as.vector(getME(nativa_interacao, "X") %*% fixef(nativa_interacao)$cond))
var_aleatoria_interacao <- VarCorr(nativa_interacao)$cond$cluster[1, 1] +
                           VarCorr(nativa_interacao)$cond$`locationName:cluster`[1, 1]

var_residual <- pi^2 / 3

tab_r2 <- data.frame(
  modelo        = c("nulo", "nativa_aditivo", "nativa_interacao"),
  var_fixos     = c(var_fixos_nulo, var_fixos_aditivo, var_fixos_interacao),
  var_aleatoria = c(var_aleatoria_nulo, var_aleatoria_aditivo, var_aleatoria_interacao),
  var_residual  = var_residual
) %>%
  mutate(R2_marginal    = var_fixos / (var_fixos + var_aleatoria + var_residual),
         R2_condicional = (var_fixos + var_aleatoria) / (var_fixos + var_aleatoria + var_residual),
         across(-modelo, ~ round(.x, 3)))
print(tab_r2, row.names = FALSE)
write.xlsx(tab_r2, here("output", "tabelas", "04c2_r2_nakagawa.xlsx"), rowNames = FALSE)

# Comparacao com o pacote performance (so para conferencia).
# Para binomial agregada (cbind), o performance usa outra variancia
# residual (aproximacao na escala observada, nao pi^2/3); por isso o R2
# condicional dele fica perto de 1. O R2 marginal fica proximo do nosso.
modelo_glmmTMB <- nativa_interacao
class(modelo_glmmTMB) <- "glmmTMB"
cat("\nperformance::r2_nakagawa (modelo nativa_interacao):\n")
print(performance::r2_nakagawa(modelo_glmmTMB))

# Leitura: o R2 condicional alto vem quase todo do efeito aleatorio de
# camera (variancia ~21; muitas cameras sem nenhum registro). O R2
# marginal e baixo: a curva de atividade e a nativa explicam pouco da
# variacao entre linhas, que e dominada pela diferenca entre cameras.


# ------------------------------------------------------------
# D. Sensibilidade
#   D1. sem P150 (97,5% de nativa, z = 3,54; 5 sucessos)
#   D2. sem P150 e P138 (os 2 pontos com muita nativa E com registro
#       de veado; 9 sucessos no total)
#   D3. nativa transformada: arco-seno da raiz da proporcao, padronizada
#       por cluster (reduz o peso dos valores extremos)
#   Em cada caso: AIC da interacao contra o nulo nos MESMOS dados
# ------------------------------------------------------------
cat("\n===== D. SENSIBILIDADE =====\n")

# D1 - sem P150
chm_sem_P150 <- chm_veado %>% filter(cluster != "P150")
nulo_sem_P150 <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                        data = chm_sem_P150)
interacao_sem_P150 <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                             data = chm_sem_P150)

# D2 - sem P150 e P138
chm_sem_P150_P138 <- chm_veado %>% filter(!cluster %in% c("P150", "P138"))
nulo_sem_P150_P138 <- chmTMB(cbind(success, failure) ~ time + (1 | cluster / locationName),
                             data = chm_sem_P150_P138)
interacao_sem_P150_P138 <- chmTMB(cbind(success, failure) ~ time * z_native_1680 + (1 | cluster / locationName),
                                  data = chm_sem_P150_P138)

# D3 - nativa transformada (arco-seno da raiz), padronizada por cluster
nativa_transf <- nativa_cluster %>%
  mutate(asin_native_1680   = asin(sqrt(native_1680 / 100)),
         z_asin_native_1680 = as.numeric(scale(asin_native_1680))) %>%
  select(cluster, z_asin_native_1680)

chm_transf <- chm_veado %>% left_join(nativa_transf, by = "cluster")
cat("Amplitude de z (original)     :", round(range(nativa_cluster$z_native_1680), 2), "\n")
cat("Amplitude de z (arco-seno)    :", round(range(nativa_transf$z_asin_native_1680), 2), "\n")

interacao_transf <- chmTMB(cbind(success, failure) ~ time * z_asin_native_1680 + (1 | cluster / locationName),
                           data = chm_transf)

cat("Convergencia:\n")
print(c(nulo_sem_P150           = nulo_sem_P150$sdr$pdHess,
        interacao_sem_P150      = interacao_sem_P150$sdr$pdHess,
        nulo_sem_P150_P138      = nulo_sem_P150_P138$sdr$pdHess,
        interacao_sem_P150_P138 = interacao_sem_P150_P138$sdr$pdHess,
        interacao_transf        = interacao_transf$sdr$pdHess))

tab_d <- data.frame(
  cenario   = c("todos os clusters",
                "sem P150",
                "sem P150 e P138",
                "nativa transformada (arco-seno)"),
  clusters  = c(n_distinct(chm_veado$cluster),
                n_distinct(chm_sem_P150$cluster),
                n_distinct(chm_sem_P150_P138$cluster),
                n_distinct(chm_transf$cluster)),
  sucessos  = c(sum(chm_veado$success),
                sum(chm_sem_P150$success),
                sum(chm_sem_P150_P138$success),
                sum(chm_transf$success)),
  AIC_nulo      = c(AIC(nulo), AIC(nulo_sem_P150), AIC(nulo_sem_P150_P138), AIC(nulo)),
  AIC_interacao = c(AIC(nativa_interacao), AIC(interacao_sem_P150),
                    AIC(interacao_sem_P150_P138), AIC(interacao_transf))
) %>%
  mutate(dAIC_interacao_nulo = round(AIC_interacao - AIC_nulo, 2),
         AIC_nulo            = round(AIC_nulo, 2),
         AIC_interacao       = round(AIC_interacao, 2))
print(tab_d, row.names = FALSE)
write.xlsx(tab_d, here("output", "tabelas", "04d_sensibilidade.xlsx"), rowNames = FALSE)

# Resultado (dAIC_interacao_nulo; negativo = interacao melhor que o nulo):
#   todos os clusters : -8,93  -> forte suporte a interacao
#   sem P150          : -2,95  -> suporte fraco (no limite de 2 AIC)
#   sem P150 e P138   : +0,32  -> o efeito DESAPARECE
#   arco-seno         : -9,16  -> a transformacao nao muda nada, porque o
#                                 P150 continua extremo (z = 3,68)
# Interpretacao: o efeito da nativa sobre o horario e sustentado por
# 9 sucessos em 2 pontos. Deve ser apresentado como resultado sugestivo,
# relatando esta analise de sensibilidade junto.


# ------------------------------------------------------------
# E. Figuras
#   Niveis de nativa: percentis 10%, 50% e 90% dos 55 clusters
#   (6%, 17% e 52%). A curva de "muita nativa" (52%) e estimada com
#   poucos pontos com registro, principalmente P150 e P138.
#
#   E1. Probabilidade de atividade (curva MARGINAL, aprox. de Zeger
#       et al. 1988). ATENCAO: o nivel absoluto nao e confiavel por causa
#       do excesso de cameras sem registro (ver script 02).
#   E2. Padrao RELATIVO: cada curva dividida pela sua soma ao longo do
#       dia (distribuicao da atividade nas 24 h, como uma densidade).
#       Mostra so a FORMA, que e o que a interacao testa.
# ------------------------------------------------------------
niveis_z   <- quantile(nativa_cluster$z_native_1680, c(0.10, 0.50, 0.90))
niveis_pct <- round(quantile(nativa_cluster$native_1680, c(0.10, 0.50, 0.90)))

var_aleatoria <- VarCorr(nativa_interacao)$cond$cluster[1, 1] +
                 VarCorr(nativa_interacao)$cond$`locationName:cluster`[1, 1]
fator_zeger <- sqrt(1 + (16 * sqrt(3) / (15 * pi))^2 * var_aleatoria)

pred <- predict(nativa_interacao, list(z_native_1680 = niveis_z)) %>%
  mutate(nivel = case_when(z_native_1680 == niveis_z[1] ~ paste0("Pouca (", niveis_pct[1], "% nativa)"),
                           z_native_1680 == niveis_z[2] ~ paste0("Media (", niveis_pct[2], "% nativa)"),
                           z_native_1680 == niveis_z[3] ~ paste0("Muita (", niveis_pct[3], "% nativa)")),
         nivel = factor(nivel, levels = c(paste0("Pouca (", niveis_pct[1], "% nativa)"),
                                          paste0("Media (", niveis_pct[2], "% nativa)"),
                                          paste0("Muita (", niveis_pct[3], "% nativa)"))),
         marginal     = plogis(fit / fator_zeger),
         marginal_lcl = plogis((fit - 1.96 * se.fit) / fator_zeger),
         marginal_ucl = plogis((fit + 1.96 * se.fit) / fator_zeger)) %>%
  mutate(relativo = marginal / sum(marginal), .by = nivel)

cores <- c("#d95f02", "#7f7f7f", "#1b9e77")

# E1 - probabilidade marginal
figura_prob <- ggplot(pred, aes(x = time, colour = nivel, fill = nivel)) +
  geom_ribbon(aes(ymin = marginal_lcl, ymax = marginal_ucl), colour = NA, alpha = 0.15) +
  geom_line(aes(y = marginal), linewidth = 1) +
  facet_wrap(~ nivel, nrow = 1) +
  scale_colour_manual(values = cores, aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Hora solar", y = "Probabilidade de atividade (por hora)",
       title = "Veado-catingueiro - vegetacao nativa em 1680 m",
       caption = "Media marginal (aprox. Zeger et al. 1988); faixa: IC 95%. Nivel absoluto pouco confiavel (excesso de zeros).") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(here("output", "figuras", "04_atividade_nativa_probabilidade.png"), figura_prob,
       width = 10, height = 4, dpi = 300)

# E2 - padrao relativo (forma da curva)
figura_rel <- ggplot(pred, aes(x = time, y = relativo, colour = nivel)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = cores) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Hora solar", y = "Proporcao da atividade diaria",
       colour = NULL,
       title = "Distribuicao da atividade ao longo do dia",
       caption = "Cada curva soma 1 ao longo das 24 h (so a forma; o nivel foi removido).") +
  theme_classic() +
  theme(legend.position = "top")
ggsave(here("output", "figuras", "04_atividade_nativa_relativa.png"), figura_rel,
       width = 7, height = 4.5, dpi = 300)

print(figura_prob)
print(figura_rel)

# Horario dos picos em cada nivel de nativa
cat("\nHorario do pico principal (hora solar):\n")
pred %>%
  slice_max(marginal, n = 1, by = nivel) %>%
  transmute(nivel, hora_pico = round(time, 1)) %>%
  print(row.names = FALSE)
