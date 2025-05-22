library('terra')
library('sf')
library("biomod2")
library("readr")
library("dplyr")
library("tidyr")
library("CoordinateCleaner")
library("usdm")
library('spThin')
library('ggplot2')
library('patchwork')
library('pROC')
library('geodata')

set.seed(112)


### setting colors to be used in the study

COLORS <- c('#caf4c6','#5dd39e','#348aa7','#525174','#513b56')

VIRI <- c('#ade8f4','#48cae4','#0096c7','#023e8a','#03045e')


### loading bioclimatic variables
nms <-
  c(
    'bio1',
    'bio2',
    'bio3',
    'bio4',
    'bio5',
    'bio6',
    'bio7',
    'bio8',
    'bio9',
    'bio10',
    'bio11',
    'bio12',
    'bio13',
    'bio14',
    'bio15',
    'bio16',
    'bio17',
    'bio18',
    'bio19'
  )

### loading bangladesh shapefile
bd <- geodata::gadm('BD', level = 0, path = getwd())

bd_sf_obj <- sf::st_as_sf(bd)

whis_data <- read_tsv('whis.csv') %>%
  select(species = "species",
         long = "decimalLongitude",
         lat = "decimalLatitude") %>%
  filter(species == 'Dendrocygna javanica')



### current bioclimatic varaibles at 2.5' resolution
world_curr <- terra::rast(
  c(
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_2.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_3.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_4.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_5.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_6.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_7.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_8.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_9.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_10.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_11.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_12.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_13.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_14.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_15.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_16.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_17.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_18.tif',
    'wc2.1_2.5m_bio/wc2.1_2.5m_bio_19.tif'
  )
)

### future bioclimatic variables using MPIESM12HR GCM at 2.5' resolution
world_2050_245 <-
  terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060.tif")
world_2050_370 <-
  terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp370_2041-2060.tif")
world_2050_585 <-
  terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2041-2060.tif")

world_2070_245 <-
  terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2061-2080.tif")
world_2070_370 <-
  terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp370_2061-2080.tif")
world_2070_585 <-
  terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2061-2080.tif")

world_2090_245 <- terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2081-2100.tif")
world_2090_370 <- terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp370_2081-2100.tif")
world_2090_585 <- terra::rast("wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2081-2100.tif")

### change names to be the same
names(world_curr) <- nms

names(world_2050_245) <- nms
names(world_2050_370) <- nms
names(world_2050_585) <- nms

names(world_2070_245) <- nms
names(world_2070_370) <- nms
names(world_2070_585) <- nms

names(world_2090_245) <- nms
names(world_2090_370) <- nms
names(world_2090_585) <- nms

### mask the current and future variables to bangladesh map
bd_curr <- terra::crop(world_curr, bd) %>% mask(bd)

bd_50_245 <- terra::crop(world_2050_245, bd) %>% mask(bd)
bd_50_370 <- terra::crop(world_2050_370, bd) %>% mask(bd)
bd_50_585 <- terra::crop(world_2050_585, bd) %>% mask(bd)

bd_70_245 <- terra::crop(world_2070_245, bd) %>% mask(bd)
bd_70_370 <- terra::crop(world_2070_370, bd) %>% mask(bd)
bd_70_585 <- terra::crop(world_2070_585, bd) %>% mask(bd)

bd_90_245 <- terra::crop(world_2090_245, bd) %>% mask(bd)
bd_90_370 <- terra::crop(world_2090_370, bd) %>% mask(bd)
bd_90_585 <- terra::crop(world_2090_585, bd) %>% mask(bd)


### reading gbif data and cleaning it

whis_data_cleaned <-
  cc_zero(whis_data, lon = "long", lat = "lat") %>% # remove zero cord
  cc_dupl(lon = "long", lat = "lat", species = "species") %>% # remove duplicate records
  cc_sea(lon = "long", lat = "lat") %>% # remove records at sea
  cc_cen(lon = 'long', lat = 'lat', species = 'species') %>% # remove records at country centroid
  cc_outl(lon = 'long', lat = 'lat', species = 'species') %>% # remove outliers
  cc_cap(lon = 'long', lat = 'lat', species = 'species') %>% # remove records at country capital
  cc_val(lon = 'long', lat = 'lat') # remove non-numeric and non-available records


### selecting target species and thinning the dataset using 5km as min thin distance
### spaitally thinning presence data to reduce sample bias
thin(
  loc.data = whis_data_cleaned,
  long.col = 'long',
  lat.col = 'lat',
  spec.col = 'species',
  thin.par = 5,
  reps = 1,
  write.files = T,
  out.dir = getwd(),
  out.base = 'whis_data_5',
  verbose = T,
  locs.thinned.list.return = F
)

sp <- read_csv('whis_data_5_thin1.csv')



### adding presence column
sp$ob <- 1

### creating spatial points data from species presence data
sp_v <-
  terra::vect(sp,
              geom = c('long', 'lat'),
              crs = crs(world_curr))

### creating background area as minimum convex polygon that includes all presence points
### and with 100km buffer around the MCP
cvb <- terra::convHull(sp_v) %>% terra::buffer(width = 100 * 1000)

### masking current climatic data to background area
curr_cvb <- terra::crop(world_curr, cvb) %>% mask(cvb)

### performing VIF analysis of 19 bioclimatic variables with 5 VIF threshold
curr_sp <- terra::extract(curr_cvb, sp_v, ID = F)
v1 <- vifstep(curr_sp, th = 5)

### excluding bioclimatic variables having more than 5 VIF score
curr_test_stack <- usdm::exclude(curr_cvb, v1)

### Formatting data and selecting two sets of pseudo absences
f_data <- BIOMOD_FormatingData(
  resp.var = as.numeric(sp$ob),
  expl.var = curr_test_stack,
  resp.xy = sp[, c("long", "lat")],
  resp.name = sp$species[1],
  filter.raster = F,
  PA.nb.rep = 2,
  PA.nb.absences = c(230, 10000),
  PA.strategy = 'disk',
  PA.dist.min = 15 * 1000,
  PA.dist.max = 150 * 1000
)

### using default options from biomod2
options <- bm_ModelingOptions(
  data.type = 'binary',
  strategy = 'default',
  models = c('RF', 'GBM', 'CTA', 'MARS', 'FDA', 'GLM', 'MAXNET')
)

### Building the model using 7 algorithms and 10 fold cross validations using each psuedoabsence datasets,
### resulting in total of (7x10) = 70 models for evaluation
mdls_eval <- BIOMOD_Modeling(
  bm.format = f_data,
  models = c('RF', 'GBM', 'CTA', 'MARS', 'FDA', 'GLM', 'MAXNET'),
  models.pa = list(
    'RF' = c('PA1'),
    'GBM' = c('PA1'),
    'CTA' = c('PA1'),
    'MARS' = c('PA1'),
    'FDA' = c('PA1'),
    'GLM' = c('PA2'),
    'MAXNET' = c('PA2')
  ),
  CV.strategy = 'kfold',
  CV.k = 10,
  CV.nb.rep = 1,
  CV.do.full.models = F,
  prevalence = 0.5,
  bm.options = options,
  var.import = 0,
  metric.eval = c('TSS', 'ROC')
)


### build final models using algorithms that have a validation AUC > 8 and TSS > 6
mdls <- BIOMOD_Modeling(
  bm.format = f_data,
  models = c('RF', 'MARS', 'MAXNET'),
  models.pa = list(
    'RF' = c('PA1'),
    'MARS' = c('PA1'),
    'MAXNET' = c('PA2')
  ),
  CV.strategy = 'kfold',
  CV.k = 2,
  CV.nb.rep = 1,
  CV.do.full.models = T,
  prevalence = 0.5,
  bm.options = options,
  var.import = 0,
  metric.eval = c('TSS', 'ROC')
)

### get model names to build the ensemble model
chosen_mdls <- BIOMOD_LoadModels(bm.out = mdls, run = 'allRun')

### build ensemble models using mean of all probability of each algos (EMmedian) and
### using 6 permutation runs to calculate variable importance
em <- BIOMOD_EnsembleModeling(
  bm.mod = mdls,
  models.chosen = chosen_mdls,
  em.by = 'all',
  em.algo = 'EMmedian',
  metric.eval = c('ROC', 'TSS'),
  var.import = 8,
  metric.select = c('ROC'),
  metric.select.thresh = 0.6,
  metric.select.dataset = 'calibration'
)

### get variable importance
varimp <- get_variables_importance(em) %>%
  select(expl.var, var.imp) %>%
  group_by(expl.var) %>%
  summarise(sd_imp = sd(var.imp), mean_imp = mean(var.imp))



### get evaluation score of the test models

mdls_test_ev_2 <- get_evaluations(mdls_eval) %>% select(algo, metric.eval, calibration, validation)


### plots ###

# 1. var importance plot

ggplot(varimp, aes(x = reorder(expl.var, mean_imp))) +
  geom_bar(mapping = aes(y = mean_imp),
           stat = 'identity',
           fill = COLORS[4]) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = 'bioclimatic variables', y = 'importance (%)') +
  theme_classic()

# 2. model evaluation plot

ggplot(mdls_test_ev_2) +
  geom_boxplot(
    mapping = aes(x = algo, y = validation, fill = metric.eval),
    color = COLORS[1]
  ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(0, 1)) +
  scale_fill_manual(values = c(COLORS[5], VIRI[4])) +
  labs(x = 'algorithm', y = 'validation score', fill = 'metric') +
  theme_minimal()


### 3. ROC curve plot

obs <- get_formal_data(em, 'resp.var')
obs[is.na(obs)] <- 0
preds <- get_predictions(em, model.as.col = T)
preds <- preds / 1000

roc_obj <- roc(obs, preds[[1]])

pROC::plot.roc(
  roc_obj,
  asp = F,
  print.auc = T,
  legacy.axes = T,
  identity.col = COLORS[2],
  print.auc.col = COLORS[4],
  col = VIRI[4]
)


### plot response curves of bio8, bio15, bio18, bio3
rc <- bm_PlotResponseCurves(bm.out = em, models.chosen = get_built_models(em))

filter(rc$tab, expl.name == 'bio3') %>%
  ggplot() +
  geom_line(mapping = aes(x = expl.val, y = pred.val),
            color = COLORS[3]) +
  scale_y_continuous(
    limits = c(0, 0.5),
    labels = c('00%', '10%', '20%', '30%', '40%', '50%')
  ) +
  labs(y = 'presence probability (%)', x = 'Isothermality') +
  theme_minimal()

filter(rc$tab, expl.name == 'bio8') %>%
  ggplot() +
  geom_line(mapping = aes(x = expl.val, y = pred.val),
            color = COLORS[3]) +
  scale_y_continuous(
    limits = c(0, 0.5),
    labels = c('00%', '10%', '20%', '30%', '40%', '50%')
  ) +
  labs(y = 'presence probability (%)', x = 'Mean Temp of Wettest Quarter (C)') +
  theme_minimal()

filter(rc$tab, expl.name == 'bio15') %>%
  ggplot() +
  geom_line(mapping = aes(x = expl.val, y = pred.val),
            color = COLORS[3]) +
  scale_y_continuous(
    limits = c(0, 0.5),
    labels = c('00%', '10%', '20%', '30%', '40%', '50%')
  ) +
  labs(y = 'presence probability (%)', x = 'Precipitation Seasonality (CV)') +
  theme_minimal()

filter(rc$tab, expl.name == 'bio18') %>%
  ggplot() +
  geom_line(mapping = aes(x = expl.val, y = pred.val),
            color = COLORS[3]) +
  scale_y_continuous(
    limits = c(0, 0.5),
    labels = c('00%', '10%', '20%', '30%', '40%', '50%')
  ) +
  labs(y = 'presence probability (%)', x = 'Precipitation of Warmest Quarter (mm)') +
  theme_minimal()





### make projections of the models on time steps 2041-2060 and 2061-2080 based on ssp245, ssp370, ssp585 scenarios
bd_curr_x <- usdm::exclude(bd_curr, v1)

bd_50_245_x <- usdm::exclude(bd_50_245, v1)
bd_50_370_x <- usdm::exclude(bd_50_370, v1)
bd_50_585_x <- usdm::exclude(bd_50_585, v1)

bd_70_245_x <- usdm::exclude(bd_70_245, v1)
bd_70_370_x <- usdm::exclude(bd_70_370, v1)
bd_70_585_x <- usdm::exclude(bd_70_585, v1)

bd_90_245_x <- usdm::exclude(bd_90_245, v1)
bd_90_370_x <- usdm::exclude(bd_90_370, v1)
bd_90_585_x <- usdm::exclude(bd_90_585, v1)

### range shift map

### plotting range change graph

ggplot(range_changes) +
  geom_line(mapping = aes(x = year, y = perc_range, colour = scenario),
            linewidth = 1.5) +
  geom_point(mapping = aes(year, y = perc_range, color = scenario),
             size = 4) +
  theme_minimal() +
  scale_color_manual(values = c('#2dc9ccff','#348aa7ff','#0b2847')) +
  scale_y_continuous(
    limits = c(0, 100),
    labels = c('00%', '25%', '50%', '75%', '100%'),
    name = 'range (%)'
  )



proj_curr <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_curr",
  new.env = bd_curr_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_50_245 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_50_245",
  new.env = bd_50_245_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_50_370 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_50_370",
  new.env = bd_50_370_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_50_585 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_50_585",
  new.env = bd_50_585_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_70_245 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_70_245",
  new.env = bd_70_245_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_70_370 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_70_370",
  new.env = bd_70_370_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_70_585 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_70_585",
  new.env = bd_70_585_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_90_245 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_90_245",
  new.env = bd_90_245_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_90_370 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_90_370",
  new.env = bd_90_370_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)

proj_90_585 <- BIOMOD_EnsembleForecasting(
  bm.em = em,
  proj.name = "whis25_f_90_585",
  new.env = bd_90_585_x,
  models.chosen = get_built_models(em),
  metric.binary = 'ROC',
  on_0_1000 = F
)



### function for making projection maps

proj_curr_map <- terra::unwrap(proj_curr@proj.out@val) %>% as.data.frame(xy = T)
colnames(proj_curr_map) <- c('x', 'y', 'presence')

proj_50_245_map <- terra::unwrap(proj_50_245@proj.out@val) %>% as.data.frame(xy = T)
colnames(proj_50_245_map) <- c('x', 'y', 'presence')

proj_50_370_map <- terra::unwrap(proj_50_370@proj.out@val) %>% as.data.frame(xy = T)
colnames(proj_50_370_map) <- c('x', 'y', 'presence')

proj_50_585_map <- terra::unwrap(proj_50_585@proj.out@val) %>% as.data.frame(xy = T)
colnames(proj_50_585_map) <- c('x', 'y', 'presence')



proj_70_245_map <- terra::unwrap(proj_70_245@proj.out@val) %>% as.data.frame(xy = T)
colnames(proj_70_245_map) <- c('x', 'y', 'presence')

proj_70_370_map <- terra::unwrap(proj_70_370@proj.out@val) %>% as.data.frame(xy = T)
colnames(proj_70_370_map) <- c('x', 'y', 'presence')

proj_70_585_map <- terra::unwrap(proj_70_585@proj.out@val) %>% as.data.frame(xy = T)
colnames(proj_70_585_map) <- c('x', 'y', 'presence')


### creating plot of maps for all the projections

gg_curr <- ggplot(data = proj_curr_map) +
  geom_raster(mapping = aes(x = x, y = y, fill = presence)) +
  scale_fill_stepsn(
    colours = VIRI,
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c('00%', '20%', '40%', '60%', '80%', '100%')
  ) +
  coord_quickmap() +
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(3, 'cm')
  )

gg_50_245 <- ggplot(data = proj_50_245_map) +
  geom_raster(mapping = aes(x = x, y = y, fill = presence)) +
  scale_fill_stepsn(
    colours = VIRI,
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c('00%', '20%', '40%', '60%', '80%', '100%'),
    guide = 'none'
  ) +
  coord_quickmap() +
  theme_void()

gg_50_370 <- ggplot(data = proj_50_370_map) +
  geom_raster(mapping = aes(x = x, y = y, fill = presence)) +
  scale_fill_stepsn(
    colours = VIRI,
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c('00%', '20%', '40%', '60%', '80%', '100%'),
    guide = 'none'
  ) +
  coord_quickmap() +
  theme_void()

gg_50_585 <- ggplot(data = proj_50_585_map) +
  geom_raster(mapping = aes(x = x, y = y, fill = presence)) +
  scale_fill_stepsn(
    colours = VIRI,
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c('00%', '20%', '40%', '60%', '80%', '100%'),
    guide = 'none'
  ) +
  coord_quickmap() +
  theme_void()

gg_70_245 <- ggplot(data = proj_70_245_map) +
  geom_raster(mapping = aes(x = x, y = y, fill = presence)) +
  scale_fill_stepsn(
    colours = VIRI,
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c('00%', '20%', '40%', '60%', '80%', '100%'),
    guide = 'none'
  ) +
  coord_quickmap() +
  theme_void()

gg_70_370 <- ggplot(data = proj_70_370_map) +
  geom_raster(mapping = aes(x = x, y = y, fill = presence)) +
  scale_fill_stepsn(
    colours = VIRI,
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c('00%', '20%', '40%', '60%', '80%', '100%'),
    guide = 'none'
  ) +
  coord_quickmap() +
  theme_void()

gg_70_585 <- ggplot(data = proj_70_585_map) +
  geom_raster(mapping = aes(x = x, y = y, fill = presence)) +
  scale_fill_stepsn(
    colours = VIRI,
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c('00%', '20%', '40%', '60%', '80%', '100%'),
    guide = 'none'
  ) +
  coord_quickmap() +
  theme_void()




### getting binary maps
b_curr <- get_predictions(proj_curr, metric.binary = 'ROC')

b_50_245 <- get_predictions(proj_50_245, metric.binary = 'ROC')
b_50_370 <- get_predictions(proj_50_370, metric.binary = 'ROC')
b_50_585 <- get_predictions(proj_50_585, metric.binary = 'ROC')

b_70_245 <- get_predictions(proj_70_245, metric.binary = 'ROC')
b_70_370 <- get_predictions(proj_70_370, metric.binary = 'ROC')
b_70_585 <- get_predictions(proj_70_585, metric.binary = 'ROC')

b_90_245 <- get_predictions(proj_90_245, metric.binary = 'ROC')
b_90_370 <- get_predictions(proj_90_370, metric.binary = 'ROC')
b_90_585 <- get_predictions(proj_90_585, metric.binary = 'ROC')



### range shift calculations of current and future potential distribution
rs_50_245 <- BIOMOD_RangeSize(b_curr, b_50_245)
rs_50_370 <- BIOMOD_RangeSize(b_curr, b_50_370)
rs_50_585 <- BIOMOD_RangeSize(b_curr, b_50_585)

rs_70_245 <- BIOMOD_RangeSize(b_curr, b_70_245)
rs_70_370 <- BIOMOD_RangeSize(b_curr, b_70_370)
rs_70_585 <- BIOMOD_RangeSize(b_curr, b_70_585)

rs_90_245 <- BIOMOD_RangeSize(b_curr, b_90_245)
rs_90_370 <- BIOMOD_RangeSize(b_curr, b_90_370)
rs_90_585 <- BIOMOD_RangeSize(b_curr, b_90_585)






range_changes <- data.frame(
  year = c(2020, 2020, 2020, 2050, 2050, 2050, 2070, 2070, 2070),
  scenario = rep(c('ssp245', 'ssp370', 'ssp585'), 3),
  range_size = c(3034, 3034, 3034, 2486, 1519, 1400, 1997, 2141, 1183)
) %>% mutate(perc_range = (range_size / 3034) * 100)


### protected area coverage

pa_calc <- function(x1) {
  x1 <- terra::as.polygons(proj_ssp245_50_mb)
  x1 <- x1[2]
  pa_x <- terra::intersect(x1, bd_pa)
  
  x_ex <- terra::expanse(x1, unit = 'km')
  pa_x_ex <- terra::expanse(pa_x, unit = 'km') %>% sum()
  
  pa_x_ex
}



### helper functions

test_disk <- function(mi, ma) {
  mi_b <-
    terra::buffer(sp_v, width = mi * 1000) %>% terra::aggregate()
  ma_b <-
    terra::buffer(sp_v, width = ma * 1000) %>% terra::aggregate()
  
  ma_i <- terra::erase(ma_b, mi_b)
  
  plot(ma_i, col = 'blue')
  plot(cvb , add = T)
}

proj_count <- function(p) {
  p_map <- terra::unwrap(p@proj.out@val) %>% as.data.frame(xy = T)
  colnames(p_map) <- c('x', 'y', 'occ')
  
  p_map <- mutate(p_map, occ = if_else(occ >= 0.6, 1, 0)) %>% count(occ)
  
  print(p_map)
}



ggplot(bd_sf_obj) +
  geom_sf(fill = 'white') +
  geom_sf(data = st_as_sf(sp_v), color = COLORS[3]) +
  coord_sf() +
  theme_void()


## f_data <- read_rds('f_data.rds')
## mdls_eval <- read_rds('mdls_eval.rds')
## mdls <- read_rds('mdls.rds')
## em <- read_rds('em.rds')
