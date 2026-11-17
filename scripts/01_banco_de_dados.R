# ============================================================
# 01 - BANCO DE DADOS PARA O chmTMB
# Veado-catingueiro (Subulo gouazoubira) - armadilhas fotograficas
#
# Entradas (data/raw/):
#   cameratrapdata.xlsx              -> cameras (esforco + paisagem)
#   cameratrapdata_native_veg.xlsx   -> usada so para corrigir o P009
#   all_species_data_solartime.xlsx  -> registros (hora solar)
#
# Saidas (data/processed/):
#   deployments.xlsx        1 linha por camera
#   observacoes_veado.xlsx  1 linha por registro de veado
#   chm_veado.xlsx / .rds   base do chmTMB (camera x hora do dia)
#
# Abra o veado_chm.Rproj antes de rodar (os caminhos usam here()).
# ============================================================


# ------------------------------------------------------------
# 0. Pacotes e funcoes do CHM
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(openxlsx)
library(here)

# Copia local do chm_helpers (Rowcliffe), com 2 correcoes:
#  - registros entre 00:00 e 01:00 eram contados como falha (R >= 4.3)
#  - type = "unimodal" nao funcionava
source(here("R", "chm_helpers.R"))


# ------------------------------------------------------------
# 1. Leitura das planilhas
# ------------------------------------------------------------
cameras_raw   <- read.xlsx(here("data", "raw", "cameratrapdata.xlsx"))
cameras_nv    <- read.xlsx(here("data", "raw", "cameratrapdata_native_veg.xlsx"))
registros_raw <- read.xlsx(here("data", "raw", "all_species_data_solartime.xlsx"))


# ------------------------------------------------------------
# 2. Correcao do P009
#   Em cameratrapdata.xlsx a retirada do P009 (15/08/2018) e ANTES da
#   instalacao, e ha registros ate 05/09/2018. A planilha native_veg
#   traz 12/09/2018 -> usamos essa data.
# ------------------------------------------------------------
retirada_nv <- cameras_nv %>%
  select(ID, camera, data_retirada_nv = data_retirada)

cameras_raw <- cameras_raw %>%
  left_join(retirada_nv, by = c("ID", "camera")) %>%
  mutate(data_retirada = if_else(toupper(ID) == "P009", data_retirada_nv, data_retirada)) %>%
  select(-data_retirada_nv)


# ------------------------------------------------------------
# 3. Tabela de cameras (deployments)
#   Hierarquia: cluster (ponto) > camera. O ID da camera inclui o
#   cluster (ex.: P015C1R), por isso e unico no estudo todo.
#
#   Datas do Excel: numero de dias desde 30/12/1899.
#   Horas do Excel: fracao do dia (0.5 = 12:00) -> x 86400 segundos.
#   Horario local de relogio; o rotulo "UTC" so evita conversao de fuso.
# ------------------------------------------------------------
deployments <- cameras_raw %>%
  transmute(
    cluster         = toupper(ID),
    camera          = camera,
    deploymentID    = paste0(cluster, camera, "R"),
    locationName    = deploymentID,
    deploymentStart = as.POSIXct(as.Date(data_instalacao, origin = "1899-12-30"), tz = "UTC") +
                        round(hora_instalacao * 86400),
    deploymentEnd   = as.POSIXct(as.Date(data_retirada, origin = "1899-12-30"), tz = "UTC") +
                        round(hora_retirada * 86400),
    dias_ativos     = as.numeric(difftime(deploymentEnd, deploymentStart, units = "days")),
    campanha        = factor(format(deploymentStart, "%Y")),
    latitude        = -as.numeric(lat),    # planilha sem o sinal negativo
    longitude       = -as.numeric(long),
    vegetacao       = vegetacao,

    # % vegetacao nativa = floresta + vegetacao herbaceo-arbustiva
    native_800  = Floresta_800  + Veg_herb_arb_800,
    native_1240 = Floresta_1240 + Veg_herb_arb_1240,
    native_1680 = Floresta_1680 + Veg_herb_arb_1680,
    native_2120 = Floresta_2120 + Veg_herb_arb_2120,
    native_2560 = Floresta_2560 + Veg_herb_arb_2560,
    native_3000 = Floresta_3000 + Veg_herb_arb_3000,

    # distancia media entre fragmentos (m)
    # NA = menos de 2 fragmentos no buffer (NAO e distancia zero)
    enn_mn_800  = enn_mn_800,
    enn_mn_1240 = enn_mn_1240,
    enn_mn_1680 = enn_mn_1680,
    enn_mn_2120 = enn_mn_2120,
    enn_mn_2560 = enn_mn_2560,
    enn_mn_3000 = enn_mn_3000,

    # NDVI medio da vegetacao nativa
    ndvi_mean_800  = ndvi_mean_800,
    ndvi_mean_1240 = ndvi_mean_1240,
    ndvi_mean_1680 = ndvi_mean_1680,
    ndvi_mean_2120 = ndvi_mean_2120,
    ndvi_mean_2560 = ndvi_mean_2560,
    ndvi_mean_3000 = ndvi_mean_3000
  )
View(deployments)

# Covariaveis padronizadas (media 0, desvio padrao 1) -> prefixo z_
# As covariaveis sao do CLUSTER (iguais nas 3 cameras), entao a
# padronizacao e feita com 1 valor por cluster.
padronizadas <- deployments %>%
  distinct(cluster, .keep_all = TRUE) %>%
  select(cluster, starts_with(c("native_", "enn_mn_", "ndvi_mean_"))) %>%
  mutate(across(-cluster, ~ as.numeric(scale(.x)), .names = "z_{.col}")) %>%
  select(cluster, starts_with("z_"))

deployments <- deployments %>%
  left_join(padronizadas, by = "cluster")


# ------------------------------------------------------------
# 4. Conferencias das cameras
# ------------------------------------------------------------
cat("\n===== CAMERAS =====\n")
cat("Cameras:", nrow(deployments), "| clusters:", n_distinct(deployments$cluster), "\n")
cat("ID de camera repetido:", sum(duplicated(deployments$deploymentID)), "\n")
cat("Fim antes do inicio:", sum(deployments$deploymentEnd <= deployments$deploymentStart), "\n")

cat("\nCameras por cluster:\n")
deployments %>% count(cluster, name = "n_cameras") %>% count(n_cameras, name = "n_clusters") %>% print()
cat("Clusters com menos de 3 cameras:\n")
deployments %>% add_count(cluster) %>% filter(n < 3) %>%
  select(cluster, locationName, campanha) %>% print()

cat("\nCovariaveis diferentes entre cameras do mesmo cluster:\n")
deployments %>%
  summarise(native = n_distinct(native_1680), enn = n_distinct(enn_mn_1680),
            ndvi = n_distinct(ndvi_mean_1680), .by = cluster) %>%
  filter(native > 1 | enn > 1 | ndvi > 1) %>% nrow() %>% print()

cat("\nDias ativos por camera:\n")
print(summary(deployments$dias_ativos))

cat("\nCameras por campanha:\n")
print(table(deployments$campanha))

cat("\nNAs de enn_mn (cameras):\n")
deployments %>% summarise(across(starts_with("enn_mn_"), ~ sum(is.na(.x)))) %>% print()


# ------------------------------------------------------------
# 5. Tabela de registros (observacoes)
#   ImageDate = data-hora de RELOGIO (serial do Excel)
#   Time      = hora SOLAR em radianos (ancorada no nascer/por do sol)
#   timestamp = data + hora solar  -> e a usada na analise
#   Se a hora solar cruzar a meia-noite em relacao ao relogio, a data
#   e ajustada (+1 ou -1 dia).
# ------------------------------------------------------------
registros <- registros_raw %>%
  transmute(
    deploymentID      = LocationName,
    especie           = Species,
    timestamp_relogio = as.POSIXct(round(ImageDate * 86400), origin = "1899-12-30", tz = "UTC"),
    hora_relogio      = as.numeric(difftime(timestamp_relogio,
                                            as.POSIXct(as.Date(timestamp_relogio), tz = "UTC"),
                                            units = "hours")),
    hora_solar        = Time * 12 / pi,
    data_solar        = case_when(hora_solar - hora_relogio < -12 ~ as.Date(timestamp_relogio) + 1,
                                  hora_solar - hora_relogio >  12 ~ as.Date(timestamp_relogio) - 1,
                                  TRUE                            ~ as.Date(timestamp_relogio)),
    timestamp         = as.POSIXct(data_solar, tz = "UTC") + round(hora_solar * 3600)
  ) %>%
  left_join(deployments %>% select(deploymentID, cluster, deploymentStart, deploymentEnd),
            by = "deploymentID") %>%
  mutate(dentro_esforco = timestamp_relogio >= deploymentStart &
                          timestamp_relogio <= deploymentEnd)

obs_veado <- registros %>%
  filter(especie == "Veado-catingueiro")


# ------------------------------------------------------------
# 6. Conferencias dos registros
# ------------------------------------------------------------
cat("\n===== REGISTROS =====\n")
cat("Registros (todas as especies):", nrow(registros), "\n")
print(table(registros$especie))
cat("Registros sem camera correspondente:", sum(is.na(registros$cluster)), "\n")
cat("Registros fora do periodo de esforco:\n")
registros %>% filter(!dentro_esforco) %>%
  select(especie, deploymentID, timestamp_relogio, deploymentStart) %>% print()
cat("Datas ajustadas (hora solar cruzou a meia-noite):",
    sum(registros$data_solar != as.Date(registros$timestamp_relogio)), "\n")

cat("\n--- Veado-catingueiro ---\n")
cat("Registros:", nrow(obs_veado),
    "| cameras:", n_distinct(obs_veado$deploymentID),
    "| clusters:", n_distinct(obs_veado$cluster), "\n")
cat("Registros fora do esforco:", sum(!obs_veado$dentro_esforco), "\n")
cat("Registros por campanha:\n")
print(table(format(obs_veado$timestamp, "%Y")))


# ------------------------------------------------------------
# 7. Base do CHM
#   make_chm_data cruza cameras e registros:
#   - divide cada dia ativo de cada camera em 24 horas (ocasioes)
#   - ocasiao com >= 1 registro = sucesso; sem registro = falha
#   - agrega por camera x hora do dia:
#       success = n de DIAS com registro naquela hora
#       failure = n de dias SEM registro naquela hora
#   covs = colunas das cameras que acompanham cada linha
# ------------------------------------------------------------
chm_veado <- make_chm_data(
  deployments,
  obs_veado,
  nBins = 24,
  covs  = c("cluster", "campanha", "vegetacao",
            names(select(deployments, starts_with(c("native_", "enn_mn_", "ndvi_mean_",
                                                    "z_native_", "z_enn_mn_", "z_ndvi_mean_")))))
)


# ------------------------------------------------------------
# 8. Conferencias da base do CHM
# ------------------------------------------------------------
# Sucessos esperados, contados direto dos registros:
# combinacoes unicas de camera x dia x hora
sucessos_esperados <- obs_veado %>%
  distinct(deploymentID, dia = as.Date(timestamp), hora = floor(hora_solar)) %>%
  nrow()

cat("\n===== BASE DO CHM =====\n")
cat("Linhas:", nrow(chm_veado), "(esperado:", nrow(deployments), "cameras x 24 =",
    nrow(deployments) * 24, ")\n")
cat("Ocasioes (success + failure):", sum(chm_veado$success + chm_veado$failure), "\n")
cat("Sucessos na base:", sum(chm_veado$success),
    "| esperados:", sucessos_esperados, "\n")
cat("Cameras com >= 1 sucesso:", n_distinct(chm_veado$locationName[chm_veado$success > 0]),
    "de", n_distinct(chm_veado$locationName), "\n")
cat("Clusters com >= 1 sucesso:", n_distinct(chm_veado$cluster[chm_veado$success > 0]),
    "de", n_distinct(chm_veado$cluster), "\n")

cat("\nSucessos por hora do dia:\n")
chm_veado %>%
  summarise(sucessos = sum(success), ocasioes = sum(success + failure), .by = timeHour) %>%
  arrange(timeHour) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nExemplo - camera P015C1R (primeiras 6 horas):\n")
chm_veado %>%
  filter(locationName == "P015C1R") %>%
  select(cluster, locationName, timeHour, timeRadian, success, failure) %>%
  head() %>%
  print()


# ------------------------------------------------------------
# 9. Salvar
# ------------------------------------------------------------
write.xlsx(deployments, here("data", "processed", "deployments.xlsx"), rowNames = FALSE)
write.xlsx(obs_veado,   here("data", "processed", "observacoes_veado.xlsx"), rowNames = FALSE)
write.xlsx(chm_veado,   here("data", "processed", "chm_veado.xlsx"), rowNames = FALSE)
saveRDS(chm_veado,     here("data", "processed", "chm_veado.rds"))

cat("\nBanco salvo em data/processed/\n")
