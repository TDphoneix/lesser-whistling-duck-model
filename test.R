library("terra")
library("biomod2")
library("readr")
library("dplyr")
library("tidyr")
library("CoordinateCleaner")
library("usdm")
library('spThin')
library('ggplot2')

setwd("/home/borshon/infinity")


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
bd <- terra::vect('aves_bd/BGD_adm/BGD_adm0.shp')

bd_pa <- terra::vect('whis/bd_pa.shp')

### current bioclimatic varaibles at 0.5' resolution
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

### future bioclimatic variables using MPIESM12HR GCM at 0.5' resolution
world_2050_245 <-
  terra::rast("worldclim_2.5/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060.tif")
world_2050_370 <-
  terra::rast("worldclim_2.5/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp370_2041-2060.tif")
world_2050_585 <-
  terra::rast("worldclim_2.5/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2041-2060.tif")

world_2070_245 <-
  terra::rast("worldclim_2.5/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2061-2080.tif")
world_2070_370 <-
  terra::rast("worldclim_2.5/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp370_2061-2080.tif")
world_2070_585 <-
  terra::rast("worldclim_2.5/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2061-2080.tif")
### change names to be the same
names(world_curr) <- nms

names(world_2050_245) <- nms
names(world_2050_370) <- nms
names(world_2050_585) <- nms

names(world_2070_245) <- nms
names(world_2070_370) <- nms
names(world_2070_585) <- nms

### mask the current and future variables to bangladesh map
bd_curr <- terra::crop(world_curr, bd) %>% mask(bd)

bd_50_245 <- terra::crop(world_2050_245, bd) %>% mask(bd)
bd_50_370 <- terra::crop(world_2050_370, bd) %>% mask(bd)
bd_50_585 <- terra::crop(world_2050_585, bd) %>% mask(bd)

bd_70_245 <- terra::crop(world_2070_245, bd) %>% mask(bd)
bd_70_370 <- terra::crop(world_2070_370, bd) %>% mask(bd)
bd_70_585 <- terra::crop(world_2070_585, bd) %>% mask(bd)

### reading gbif data and cleaning it
whis_bd <- read_tsv("whis_raw.csv") %>%
  select(species = "species",
         long = "decimalLongitude",
         lat = "decimalLatitude")

whis_bd_cleaned <-
  cc_zero(whis_bd, lon = "long", lat = "lat") %>% # remove zero cord
  cc_dupl(lon = "long", lat = "lat", species = "species") %>% # remove duplicate records
  cc_sea(lon = "long", lat = "lat") %>% # remove records at sea
  cc_cen(lon = 'long', lat = 'lat', species = 'species') %>% # remove records at country centroid
  cc_outl(lon = 'long', lat = 'lat', species = 'species') %>% # remove outliers
  cc_cap(lon = 'long', lat = 'lat', species = 'species') %>% # remove records at country capital
  cc_val(lon = 'long', lat = 'lat') # remove non-numeric and non-available records


### selecting target species and thinning the dataset using 5km as min thin distance
### spaitally thinning presence data to reduce sample bias
thin(
  loc.data = whis_bd_cleaned,
  long.col = 'long',
  lat.col = 'lat',
  spec.col = 'species',
  thin.par = 5,
  reps = 10,
  write.files = T,
  out.dir = getwd(),
  out.base = 'w_sp_5',
  verbose = T,
  locs.thinned.list.return = F
)

sp <- read_csv('w_sp_5_thin1.csv')



### adding presence column
sp$ob <- 1

### creating spatial points data from species presence data
sp_v <-
  terra::vect(sp,
              geom = c('long', 'lat'),
              crs = crs(world_curr))

### creating background area using 100km buffer around each presence points
cvb <- terra::convHull(sp_v) %>% terra::buffer(width = 100 * 1000)

### masking current climatic data to background area
curr <- terra::crop(world_curr, cvb) %>% mask(cvb)

### performing VIF analysis of 19 bioclimatic variables with 5 VIF threshold
curr_sp <- terra::extract(curr, sp_v, ID = F)
v1 <- vifstep(curr_sp, th = 5)

### excluding bioclimatic variables having more than 5 VIF score
curr_test_stack <- usdm::exclude(curr, v1)

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
options <- BIOMOD_ModelingOptions()

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


### evaluate the models based on cross validation runs (boxplot and cross plot)
ev <- get_evaluations(mdls_eval)

eval <- bm_PlotEvalBoxplot(
  mdls_eval,
  dataset = 'validation',
  group.by = c('algo', 'algo'),
  do.plot = F
)

eval$plot + expand_limits(x = c(0, 1), y = c(0, 1))

evalm <- bm_PlotEvalMean(
  bm.out = mdls_eval,
  metric.eval = c('ROC', 'TSS'),
  dataset = 'validation',
  group.by = 'algo',
  do.plot = F
)

evalm$plot + expand_limits(x = c(0, 1), y = c(0, 1))


### build final models using algorithms that have a validation AUC > 8 and TSS > 6
mdls <- BIOMOD_Modeling(
  bm.format = f_data,
  models = c('RF', 'GBM', 'CTA', 'MARS', 'MAXNET'),
  models.pa = list(
    'RF' = c('PA1'),
    'GBM' = c('PA1'),
    'CTA' = c('PA1'),
    'MARS' = c('PA1'),
    'MAXNET' = c('PA2')
  ),
  CV.strategy = 'kfold',
  CV.k = 10,
  CV.nb.rep = 1,
  CV.do.full.models = T,
  prevalence = 0.5,
  bm.options = options,
  var.import = 0,
  metric.eval = c('TSS', 'ROC')
)

### get model names to build the ensemble model
chosen_mdls <- BIOMOD_LoadModels(bm.out = mdls,
                                 run = 'allRun')

### build ensemble models using mean of all probability of each algos (EMmedian) and
### using 6 permutation runs to calculate variable importance
em <- BIOMOD_EnsembleModeling(
  bm.mod = mdls,
  models.chosen = chosen_mdls,
  em.by = 'all',
  em.algo = 'EMmedian',
  metric.eval = c('ROC', 'TSS'),
  var.import = 6,
  metric.select = c('ROC'),
  metric.select.thresh = 0.6,
  metric.select.dataset = 'calibration'
)


### get variable importance
varimp <- get_variables_importance(em) %>%
  select(expl.var, var.imp) %>%
  group_by(expl.var) %>%
  summarize(vimp = mean(var.imp))

### plot response curves of each variable
rc <- bm_PlotResponseCurves(bm.out = em,
                            models.chosen = get_built_models(em))


### make projections of the models on time steps 2041-2060 and 2061-2080 based on ssp245, ssp370, ssp585 scenarios
bd_curr_x <- usdm::exclude(bd_curr, v1)

bd_50_245_x <- usdm::exclude(bd_50_245, v1)
bd_50_370_x <- usdm::exclude(bd_50_370, v1)
bd_50_585_x <- usdm::exclude(bd_50_585, v1)

bd_70_245_x <- usdm::exclude(bd_70_245, v1)
bd_70_370_x <- usdm::exclude(bd_70_370, v1)
bd_70_585_x <- usdm::exclude(bd_70_585, v1)


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

### getting binary maps
b_curr <- get_predictions(proj_curr, metric.binary = 'ROC')

b_50_245 <- get_predictions(proj_50_245, metric.binary = 'ROC')
b_50_370 <- get_predictions(proj_50_370, metric.binary = 'ROC')
b_50_585 <- get_predictions(proj_50_585, metric.binary = 'ROC')

b_70_245 <- get_predictions(proj_70_245, metric.binary = 'ROC')
b_70_370 <- get_predictions(proj_70_370, metric.binary = 'ROC')
b_70_585 <- get_predictions(proj_70_585, metric.binary = 'ROC')


### range shift calculations of current and future potential distribution
rs_50_245 <- BIOMOD_RangeSize(b_curr, b_50_245)
rs_50_370 <- BIOMOD_RangeSize(b_curr, b_50_370)
rs_50_585 <- BIOMOD_RangeSize(b_curr, b_50_585)

rs_70_245 <- BIOMOD_RangeSize(b_curr, b_70_245)
rs_70_370 <- BIOMOD_RangeSize(b_curr, b_70_370)
rs_70_585 <- BIOMOD_RangeSize(b_curr, b_70_585)


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