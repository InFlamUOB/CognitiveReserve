#try and get libraries working  - just 4.2.0 and restart session!!

#SessionInfoa()
#getOption("defaultPackages")
#.libPaths()


##!/usr/bin/env Rscript
#.libPaths( c( .libPaths(), "/rds/homes/l/lxb732/R/library/4.0.0/haswell") )

#.libPaths()
#.libPaths( c( .libPaths(), "/rds/homes/b/bravol/R/library/4.2.0/EL8-ice") )
#.libPaths()
#
#myPaths <- .libPaths()  
#myPaths <- c(myPaths[3], myPaths[1] )  # switch them
#.libPaths(myPaths)


#setwd("/rds/projects/g/gkoutosg-variant-prediction/Laura/Magda")
#setwd("/castles/nr/projects/c/chechlmy-brain-ageing-ukbiobank/Peter")

library(dplyr)
library(workflowsets)
library(tidymodels)
library(data.table)
library(ggplot2)
library(finetune)
library(parsnip)
library(baguette)
library(DALEX)
library(fastshap)
library(vip)
library(NeuralNetTools)
library(stringr)
library(patchwork)

options(bitmapType='cairo')


print("here")

#load("/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/Magda/FinalPipeline_20211024_2033/ModelingWorkflows.RData")
#load("/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/Magda/FinalPipeline_20211024_2033/Modeling2data3mods.RData")


mainDir <-"/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/Magda"
NameRun <- 'Oct' #HackathonLong before name
subDir <- paste0(sub('\\..*', '', NameRun), format(Sys.time(), '_%Y%m%d_%H%M'))
dir.create(file.path(mainDir, subDir))
IncludeFigHere <- file.path(mainDir, subDir)

############################# 

WithCogAndSNPs4ModeAPOE <- read.csv("/rds/projects/c/chechlmy-brain-ageing-ukbiobank/Laura/WithCogAndSNPs4ModeAPOE.csv")

############################ mice

# ISARIC WHO CCP-UK study: 4C Mortality Score
# Multiple imputation
# 03_mice.R
# Centre for Medical Informatics, Usher Institute, University of Edinburgh 2020

# 1. Missing data inspection, description and characterisation 
# 2. Basic MICE function

# --------------------------------------------------------------------

# Variables for mice (Multivariate Imputation by Chained Equations)
library(finalfit)
library(dplyr)
library(mice)
library(tidylog)



# Missing data inspection
WithCogAndSNPs4ModeAPOE %>%
  missing_plot(.)

############################# Visualize Features - decide which combiations of efatures we go for, aswell as main outcome. 


#IncludeFigHere <- "/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/Magda/August_20230822_1743/"
Model1aPivot <- WithCogAndSNPs4ModeAPOE %>%
  dplyr::select(-Age2) %>%
  pivot_longer(-c(FID), names_to = "Features", values_to = "Value") %>%
  mutate(Type = case_when(startsWith(Features, "WM") | startsWith(Features, "c") ~ "Networks", 
                          startsWith(Features, "rs") | startsWith(Features, "e") ~ "Genes", 
                          Features %in% c("ClinicalAge","Clinical_Sex")  ~ "Demographics", 
                          Features == "Clinical_EducationMeasurement" ~ "Education",
                          TRUE ~ "Cognition"))

library(RColorBrewer)
myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(as.factor(Model1aPivot$Type))
colScale <- scale_fill_manual(name = "Type",values = myColors)

Networks <- Model1aPivot %>%
  filter(Type == "Networks" | Features %in% c("ClinicalAge", "GeneralCog", "Clinical_EducationMeasurement") ) # %>% 
# mutate(Type2 = str_replace(Type, "Networks", "Brain Imaging"))
# mutate(Type = ifelse(Type == "Networks", "Brain Imaging", .)) %>%
# mutate(Type = ifelse(Type == "Networks", "Brain Imaging", .)) %>%
#   mutate(Type = ifelse(Type == "Networks", "Brain Imaging", .)) %>%  

#  General cognitive ability score (GCAS) 


NetworksPlot <- Networks %>%
  ggplot(aes(Value)) +
  facet_wrap(~Features, scales = "free") +
  geom_histogram(aes(fill = Type)) + 
  theme_bw() +  colScale


pdf(paste0(IncludeFigHere, "/ContinousPLot.pdf"), 16,  9)
print(NetworksPlot)
dev.off()

Genes <- Model1aPivot %>%
  filter( Features %in% c( "Clinical_Sex") | Type == "Genes") %>%
  mutate( Value = as.character(Value)) %>%
  count( Features, Value ,sort = TRUE) %>%
  group_by(Features) %>%
  mutate(Sum = sum(n), 
         Freq = round(n/Sum, 2 )) %>%
  mutate(Type = case_when(Features %in% c("ClinicalAge","Clinical_Sex") ~ "Demographics", TRUE ~ "Genes")) %>%
  
  ungroup() %>%
  ggplot(aes(Freq, Value )) +
  geom_col(aes(fill = Type)) +
  facet_wrap(~Features, scales = "free_y") + 
  theme_bw() +   colScale

pdf(paste0(IncludeFigHere, "/FactorPLot.pdf"), 13,  9)
print(Genes)
dev.off()

pdf(paste0(IncludeFigHere, "/AllPlots.pdf"), 13,  15)
NetworksPlot/Genes
dev.off()

######################################  Preprocessing


data <- WithCogAndSNPs4ModeAPOE %>% 
  dplyr::select(-c(FID, Age2)) %>%
  dplyr::mutate_at(vars(Clinical_Sex, rs80136977_A:e2e3), as.factor)

col_y <- "GeneralCog" #outcome
plotType <- "D" #plot format
y <- data[[col_y]]

.seed <- 11
set.seed(.seed)

splits <- initial_split(data, strata = col_y, prop = 5/6)
data_trn <- splits %>% rsample::training()
data_tst <- splits %>% rsample::testing()

######################################  Model Building 1 (Recipes)

linear_reg_spec <- 
  linear_reg() %>% 
  set_engine("lm")

lasso_reg_spec <- 
  linear_reg(penalty = tune(), mixture = 1) %>% 
  set_engine("glmnet")

nnet_spec <- 
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>% 
  set_engine("nnet", MaxNWts = 2600 ) %>% #(MaxNWts=84581)
  set_mode("regression")

rf_spec <- 
  rand_forest( trees = 1000) %>% #mtry = tune(), min_n = tune(), 
  set_engine("ranger", importance = 'permutation') %>% 
  set_mode("regression")

xgb_spec <- 
  boost_tree(tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(), 
             min_n = tune(), sample_size = tune(), trees = tune()) %>% #tree_depth = tune(),#learn_rate = tune(), loss_reduction = tune(),
  set_engine("xgboost") %>% 
  set_mode("regression")


############################ Model Building 2 (Recipes)
# 14 different models (combinations of each)


AllRec <- recipe(GeneralCog ~ ., data = data_trn)
PreparedPreProc <- AllRec %>% prep()
AllRecJ <- juice(PreparedPreProc) #training data - basis of the rest

AllRec2 <- AllRec %>%
  step_normalize(all_numeric()) %>% #had skip = TRUE : took it out
  step_dummy(all_nominal(),  one_hot = FALSE)  %>% #change all from true to false
  step_zv(all_predictors())
PreparedPreProc <- AllRec2 %>% prep()
AllRec2J <- juice(PreparedPreProc)  # All variables dummified

#----------

ImagingRec <- AllRec %>% 
  step_rm(starts_with("Clinical"),starts_with("e"), starts_with("rs")) %>%
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- ImagingRec %>% prep()
Imaging <- juice(PreparedPreProc) #Imaging dataset ( only imaging information and outcome)

ImagingGenesRec <- AllRec %>% 
  step_rm(starts_with("Clinical")) %>%
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- ImagingGenesRec %>% prep()
ImagingGenes <- juice(PreparedPreProc) #Imaging dataset ( only imaging information and outcome)

ImagingRecEducation <- AllRec %>% 
  step_rm("Clinical_Sex", "ClinicalAge", starts_with("e"), starts_with("rs")) %>%
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- ImagingRecEducation %>% prep()
ImagingEducation <- juice(PreparedPreProc) #Imaging dataset ( only imaging information and outcome)

ImagingRecClinicalEducation <- AllRec %>% 
  step_rm(starts_with("e"), starts_with("rs")) %>%
  step_normalize(all_numeric()) %>% # skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- ImagingRecClinicalEducation %>% prep()
ImagingClinicalEducation <- juice(PreparedPreProc) #Imaging dataset ( only imaging information and outcome)


#before no clinical x imaging? CHECK!!!!!!!!!!

ImagingRecClinical <- AllRec %>% 
  step_rm(starts_with("c", ignore.case = FALSE),starts_with("e"), starts_with("rs"),"Clinical_EducationMeasurement") %>% 
  step_normalize(all_numeric()) %>% # skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- ImagingRecClinicalEducation %>% prep()
ImagingClinicalEducation <- juice(PreparedPreProc) #Imaging dataset ( only imaging information and outcome)


ImagingRecClinicalGenes <- AllRec %>% 
  step_rm("Clinical_EducationMeasurement") %>%
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- ImagingRecClinicalGenes %>% prep()
ImagingClinicalGenes <- juice(PreparedPreProc) #Imaging dataset ( only imaging information and outcome)

ImagingRecEducationGenes <- AllRec %>% 
  step_rm("Clinical_Sex", "ClinicalAge") %>%
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- ImagingRecEducationGenes %>% prep()
ImagingEducationGenes <- juice(PreparedPreProc) #Imaging dataset ( only imaging information and outcome)

#----------

GenesRec <- AllRec %>%   
  step_rm(starts_with("c"),starts_with("WM")) %>% 
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE)  %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- GenesRec %>% prep()
Genes <- juice(PreparedPreProc) #Genes dataset ( only genes information and outcome)

GenesRecEducation <- AllRec %>%   
  step_rm(starts_with("c",ignore.case = FALSE), starts_with("WM"), "Clinical_Sex", "ClinicalAge" ) %>% # no e (APOE)
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE)  %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- GenesRecEducation %>% prep()
GenesEducation <- juice(PreparedPreProc) #Genes dataset ( only genes information and outcome)


GenesRecEducationClinical <- AllRec %>%   
  step_rm(starts_with("c", ignore.case = FALSE),starts_with("WM")) %>% 
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE)  %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- GenesRecEducationClinical %>% prep()
GenesEducationClinical <- juice(PreparedPreProc) #Genes dataset ( only genes information and outcome)

GenesRecClinical <- AllRec %>%   
  step_rm(starts_with("c", ignore.case = FALSE),starts_with("WM"),"Clinical_EducationMeasurement") %>% 
  step_normalize(all_numeric()) %>% # skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE)  %>%
  step_zv(all_predictors()) #added
PreparedPreProc <- GenesRecClinical %>% prep()
GenesClinical <- juice(PreparedPreProc) #Genes dataset ( only genes information and outcome)

#----------

ClinicalRec <- AllRec %>%  
  step_rm(starts_with("c",ignore.case = FALSE ), starts_with("WM"),"Clinical_EducationMeasurement", starts_with("rs"), starts_with("e")) %>%
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) 
PreparedPreProc <- ClinicalRec %>% prep()
Clinical <- juice(PreparedPreProc)  #Just clinical information (age/sex..) 

EducationAndClinicalRec <- AllRec %>% 
  step_rm(starts_with("c",ignore.case = FALSE ), starts_with("WM"), starts_with("rs"), starts_with("e")) %>%
  step_normalize(all_numeric()) %>% #skip = TRUE
  step_dummy(all_nominal(),  one_hot = FALSE) #in theory xgboost does not accept this very well?
PreparedPreProc <- EducationAndClinicalRec %>% prep()
EducationClinical <- juice(PreparedPreProc) #Just clinical information (age/sex..) and education 

#----------

EducationRec <- AllRec %>% 
  step_rm(starts_with("c",ignore.case = FALSE ), "ClinicalAge", "Clinical_Sex", starts_with("e"),starts_with("WM"), starts_with("rs") ) %>%
  step_normalize(all_numeric())  #skip = TRUE
PreparedPreProc <- EducationRec %>% prep()
Education <- juice(PreparedPreProc) #Just clinical information (age/sex..) and education 


# before we had dummy before normalize and so all was numeric
################


recipe_list <- 
  list( AllRec2J = AllRec2,
        Imaging = ImagingRec, ImagingEducationGenes = ImagingRecEducationGenes, 
        ImagingClinicalGenes = ImagingRecClinicalGenes, ImagingEducation = ImagingRecEducation,
        ImagingClinicalEducation = ImagingRecClinicalEducation, ImagingGenes = ImagingGenesRec,
        EducationClinical = EducationAndClinicalRec,  Education = EducationRec,
        Clinical = ClinicalRec, 
        GenesEducation = GenesRecEducation, 
        Genes = GenesRec, 
        GenesEducationClinical = GenesRecEducationClinical,
        GenesClinical = GenesRecClinical
  ) 

#recipe_list <- 
#  list( AllRec2J = AllRec2)
model_list <- 
  list( nnet = nnet_spec, lasso = lasso_reg_spec , ranger = rf_spec, glm = linear_reg_spec,xgboost = xgb_spec )# xgboost = xgb_spec not working properly

#model_list <- 
#  list( lasso = lasso_reg_spec, nnet = nnet_spec, glm = linear_reg_spec)# glm not working!!- need newsx?

#model_list <- 
#  list( lasso = lasso_reg_spec,glm = linear_reg_spec)# glm not working!!- need newsx?

#########################


model_set <- workflow_set(preproc = recipe_list, models = model_list, cross = T)

model_set <- model_set %>% 
  anti_join(tibble(wflow_id = c( "Education_lasso", # cannot build a single predictor model 
                                 by = "wflow_id"))) #"AllInteractionsGenesImagJ_nnet", "AllInteractionsJ_nnet",

Run <- "Trial"
Total <- paste0(Run,"_", Sys.time())


#train_resamples <- bootstraps(uni_train, times = 5) #cannot reduce bootstrap to 2!!!
cv_splits <- rsample::vfold_cv(data_trn, v = 100, strata = col_y)

library(doParallel)
#registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
registerDoParallel(cores=60) #20


multi_metric2 <- metric_set(rmse, rsq, mae)

race_ctrl <-
  control_race(
    save_pred = TRUE,
    parallel_over = "everything",
    save_workflow = TRUE
  )

race_results <-
  model_set %>%
  workflow_map(
    "tune_race_anova",
    seed = 1503,
    resamples = cv_splits,
    grid = 10,
    control = race_ctrl, 
    verbose = TRUE,
    metrics = multi_metric2
  )


#saveRDS(race_results,paste0(Total, "race.RDS"))
save(race_results, model_set,IncludeFigHere,AllRec2J,  file = paste0(IncludeFigHere,"/Modeling2data3mods.RData"))

#load("/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/Magda/Oct_20231026_1538/Modeling2data3mods.RData")

get_model <- function(x) {
  pull_workflow_fit(x) %>% tidy(conf.int = TRUE) #tidy conf int not working? 
}
ctrl <- control_resamples(extract = get_model, save_pred = TRUE)

bootstraps <- 100 #100
Boots <-  bootstraps(data_trn, times = bootstraps)

############ Feature importance


choose_f_predict <- function(i) {
  #engine <- match.arg(engine)
  
  #engine<-  str_split(i, "_")[[1]][[2]]
  
  engine <- i
  
  f_generic <- function(object, newdata) as.numeric(predict(object, newdata = newdata))
  fs <-
    list(
      'xgboost' = f_generic,
      'ranger' = function(object, newdata) as.numeric(predict(object, data = newdata)$predictions),
      'glm' = f_generic,
      'nnet' = f_generic,
      # Choosing no penalty.
      'lasso' = function(object, newdata) as.numeric(predict(object, newx = newdata, s = 0))
    )
  fs[[engine]]
}


choose_data_setX <- function(i) {
  
  data_trn_jui <- eval(parse(text = str_split(i, "_")[[1]][[1]]))
  
  x_trn_jui <-  data_trn_jui[, setdiff(names(data_trn_jui), col_y)] %>% as.matrix()
  
  
}

choose_data_setY <- function(i) {
  
  data_trn_jui <-  eval(parse(text = str_split(i, "_")[[1]][[1]]))
  
  y_trn_jui <- data_trn_jui[[col_y]]
  
}

best_results <- list()
boosting_test_resultsWF <- list()
boosting_test_results <- list()
vi_rnks <- list()
hh <- list()

VarImp <- race_results$wflow_id %>% 
  setdiff(c("Education_nnet", "Education_glm", "Education_ranger", "Education_xgboost", "Education_lasso")) ####add all education related results here, why was lasso not here?

VarImp <- VarImp[!grepl("glm", VarImp)]

set.seed(.seed)

#library(doMC)
library(doRNG)

registerDoRNG()
#seed <- doRNGseed()
#result <- foreach(i = race_results$wflow_id) %dopar% { #, .combine = cbind
#result <- foreach(i = VarImp, .combine = cbind) %dopar% { #, .combine = cbind
result <- foreach(i = VarImp, .options.RNG=1234, .combine = cbind) %dorng% { #, .combine = cbind
  
 # doRNGseed(seed)
  print(paste0("....................", i))
  
  best_results[[i]] <- 
    race_results %>% 
    extract_workflow_set_result(i) %>% 
    select_best(metric = "rmse")
  
  boosting_test_resultsWF[[i]] <- 
    race_results %>% 
    extract_workflow(i) %>%
    finalize_workflow(best_results[[i]]) %>% #ad
    fit(data_trn) %>% 
    workflows::pull_workflow_fit()
  
  
  # boosting_test_results[[i]] <- race_results %>% 
  #   extract_workflow(i) %>%
  #   finalize_workflow(best_results[[i]]) %>% #add here the boostrapping? 
  #   fit_resamples(resamples = Boots, control = ctrl)
  #
  
  
  
  
  ######################
  #####################
  
  
  res <-
    vip::vip(
      method = 'model',
      object = boosting_test_resultsWF[[i]]$fit, 
      num_features = choose_data_setX(i) %>% ncol()
    ) %>% 
    pluck('data') %>% 
    # Will get a "Sign" column when using the default `method = 'model'`.
    rename(var = Variable, imp = Importance)
  
  if(any(names(res) == 'Sign')) {
    res <-
      res %>% 
      mutate(dir = ifelse(Sign == 'POS', +1L, -1L)) %>% 
      mutate(imp = dir * imp)
  }
  
  
  vi_vip_model <- res
  
  ####### PERMUTE
  
  f_predict <- choose_f_predict(str_split(i, "_")[[1]][[2]])
  
  resP <-
    vip::vip(
      method = 'permute',
      object = boosting_test_resultsWF[[i]]$fit, 
      num_features = choose_data_setX(i) %>% ncol(),
      metric = "rmse", 
      pred_wrapper = f_predict, #works for glm and glmnet 
      train = choose_data_setX(i) ,
      target = choose_data_setY(i)
    ) %>% 
    pluck('data') %>% 
    # Will get a "Sign" column when using the default `method = 'model'`.
    rename(var = Variable, imp = Importance)
  
  if(any(names(resP) == 'Sign')) {
    resP <-
      resP %>% 
      mutate(dir = ifelse(Sign == 'POS', +1L, -1L)) %>% 
      mutate(imp = dir * imp)
  }
  
  vi_vip_permute <- resP
  
  ####### SHAP
  
  
  resP <-
    vip::vip(
      method = 'shap',
      object = boosting_test_resultsWF[[i]]$fit, 
      num_features = choose_data_setX(i) %>% ncol(),
      train = choose_data_setX(i), 
      pred_wrapper = f_predict, 
      #exact = TRUE - only xgboost and lm
    ) %>% 
    pluck('data') %>% 
    # Will get a "Sign" column when using the default `method = 'model'`.
    rename(var = Variable, imp = Importance)
  
  if(any(names(resP) == 'Sign')) {
    resP <-
      resP %>% 
      mutate(dir = ifelse(Sign == 'POS', +1L, -1L)) %>% 
      mutate(imp = dir * imp)
  }
  
  vi_vip_shap <- resP
  
  
  ####### DALEX
  
  if (  engine<-  str_split(i, "_")[[1]][[2]] != "xgboost"){
    
    expl_dalex <- 
      DALEX::explain(
        boosting_test_resultsWF[[i]]$fit,  
        data = as.data.frame(choose_data_setX(i)),
        y = choose_data_setY(i), 
        verbose = FALSE
      )
    
  } else {
    
    expl_dalex <- 
      DALEX::explain(
        boosting_test_resultsWF[[i]]$fit,  
        data = (choose_data_setX(i)),   #xgboost cannot have the dataframe
        y = choose_data_setY(i), 
        verbose = FALSE
      )
    
  }
  
  vi_dalex_init <- 
    expl_dalex %>% 
    DALEX::variable_importance(
      type = 'difference',
      loss_function = DALEX::loss_root_mean_square, 
      n_sample = NULL
    )
  vi_dalex_init
  
  vi_dalex <-
    vi_dalex_init %>% 
    as_tibble() %>% 
    filter(permutation == 0) %>% 
    mutate(
      imp = abs(dropout_loss) / max(abs(dropout_loss))
    ) %>% 
    select(var = variable, imp) %>%
    filter(!(var %in% c('_baseline_', '_full_model_'))) %>% 
    arrange(desc(imp))
  
  
  vi_rnks[[i]] <-
    list(
      vip_model = vi_vip_model,
      vip_permute = vi_vip_permute,
      vip_shap = vi_vip_shap,
      # fastshap = vi_fastshap,
      dalex = vi_dalex
    ) %>% 
    map_dfr(bind_rows, .id = 'src') %>% 
    group_by(src) %>% 
    mutate(imp_abs = abs(imp)) %>% 
    mutate(imp_abs_norm = imp_abs / sum(imp_abs)) %>% 
    select(var, imp, imp_abs, imp_abs_norm) %>% 
    mutate(rnk = row_number(desc(imp_abs))) %>% 
    ungroup() %>%
    add_column(Model = str_split(i, "_")[[1]][[2]]) %>%
    add_column(Dataset = str_split(i, "_")[[1]][[1]])
  
  #hh[[i]] <- list(best_results[[i]],boosting_test_resultsWF[[i]],boosting_test_results[[i]],  vi_rnks[[i]])
  #if (i == "GenesClinical_ranger"){
  #  save(vi_rnks[[i]], file = paste0(IncludeFigHere,"/vi_ranksCheck_",i,".RData"))
 # }
  
  vi_rnks
  
}

# have to include engine here? 

print("Finished big loop ............... now saving ")

#save(vi_rnks, model_set,  file = paste0(IncludeFigHere,"/ModelingWorkflows.RData"))
save(result, model_set,  file = paste0(IncludeFigHere,"/ModelingWorkflows.RData"))

prettify_engine_col <- function(data) {
  res <- data %>% mutate_at(vars(Model), ~sprintf('{%s}', Model))
}

factor_src <- function(x) {
  ordered(x, levels = c('vip_model', 'vip_shap', 'vip_permute', 'dalex'))
}

names(result) <- VarImp

vi_rnks2 <- result %>%
  bind_rows() %>%
  as.data.frame()

Threshold <- 10

#IncludeFigHere <- "/rds/projects/g/gkoutosg-variant-prediction/Laura/Magda/FinalPipeline_20211024_2033"


for (j in unique(vi_rnks2$Dataset)){
  
  viz <-
    vi_rnks2 %>% 
    filter(Dataset == j) %>%
    group_by(var) %>% 
    mutate(rnk_mean = rnk %>% mean(na.rm = TRUE)) %>% 
    ungroup() %>% 
    #mutate_at(vars(var), ~forcats::fct_reorder(., -rnk_mean)) %>% 
    #ungroup() %>% 
    #prettify_engine_col() %>% 
    #mutate_at(vars(src), ~ordered(., levels = c('vip_model', 'vip_shap', 'vip_permute', 'dalex'))) %>% 
    mutate(lab = sprintf('%2d (%s)', rnk, scales::percent(imp_abs_norm, accuracy = 1, width = 2, justify = 'right'))) %>% 
    filter(rnk < Threshold) %>%
    ggplot() +
    aes(x = src, y = var) +
    geom_tile(aes(fill = rnk), alpha = 0.5, show.legend = F) +
    geom_text(aes(label = lab)) +
    scale_fill_viridis_c(direction = -1, option = "D", na.value = 'white') +
    theme_minimal(base_family = '') +
    facet_wrap(~Model) +
    theme(
      plot.title.position = 'plot',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = 'bold'),
      plot.subtitle = ggtext::element_markdown(),
    ) +
    labs(x = NULL, y = NULL)
  
  pdf(paste0(IncludeFigHere,"/NewHeatmaps_",j,".pdf"), 15, 12)
  
  print(viz)
  
  dev.off()
  
  ggsave(paste0(IncludeFigHere,"/NewHeatmapSave_",j,".pdf"), plot = viz, device = "pdf")
  #ggsave(paste0(IncludeFigHere,"/HeatmapSave_",j,".png"), plot = viz, device = "png")
  
}


Pred <- collect_predictions(race_results) 

#Likelihood ratio test? 

Perf <- collect_metrics(race_results) %>% 
  separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") %>%
  filter(.metric == "rmse") %>% 
  group_by(wflow_id) %>% 
  filter(mean == max(mean)) %>% 
  group_by(model) %>% 
  dplyr::select(-.config) %>% 
  distinct() %>%
  ungroup() %>% 
  mutate(Workflow_Rank =  row_number(mean),
         .metric = str_to_upper(.metric)) %>%
  mutate(ContainsClinical = ifelse(str_detect(Recipe, "Clinical"), 1,0)) %>%
  mutate(ContainsImaging = ifelse(str_detect(Recipe, "Imaging"), 1,0)) %>%
  mutate(ContainsGenes = ifelse(str_detect(Recipe, "Genes"),1,0)) %>%
  filter(Model_Type != "xgboost")


PlotClinical <- Perf %>%
  filter(ContainsClinical == 1 | Recipe == "AllRec2J" ) %>%
  mutate(Recipe = str_remove(Recipe, "Clinical")) %>%
  ggplot(aes(x=Workflow_Rank, y = mean, shape = Model_Type, color = Recipe)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-std_err, ymax = mean+std_err)) +
  theme_minimal()+
  #scale_colour_manual(values = cols) +
  labs(title = "Performance Comparison of Models - CLINICAL", x = "Model Rank", y = "RMSE", color = "Data", shape = "Model") #+
#scale_x_continuous(breaks = seq(from = 1, to = 30, by = 2))

pdf(paste0(IncludeFigHere,"/PlotClinical.pdf"), 6,5 )
print(PlotClinical)
dev.off()


PlotImaging <- Perf %>%
  filter(ContainsImaging == 1 | Recipe == "AllRec2J" ) %>%
  mutate(Recipe = str_remove(Recipe, "Imaging")) %>%
  ggplot(aes(x=Workflow_Rank, y = mean, shape = Model_Type, color = Recipe)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-std_err, ymax = mean+std_err)) +
  theme_minimal()+
  #scale_colour_manual(values = cols) +
  labs(title = "Performance Comparison of Models - Imaging", x = "Model Rank", y = "RMSE", color = "Data", shape = "Model") #+
#scale_x_continuous(breaks = seq(from = 1, to = 30, by = 2))


pdf(paste0(IncludeFigHere,"/PlotImaging.pdf"), 6,5 )
print(PlotImaging)
dev.off()


PlotGenes <- Perf %>%
  filter(ContainsGenes == 1 | Recipe == "AllRec2J" ) %>%
  mutate(Recipe = str_remove(Recipe, "Genes")) %>%
  ggplot(aes(x=Workflow_Rank, y = mean, shape = Model_Type, color = Recipe)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-std_err, ymax = mean+std_err)) +
  theme_minimal()+
  #scale_colour_manual(values = cols) +
  labs(title = "Performance Comparison of Models - Genes", x = "Model Rank", y = "RMSE", color = "Data", shape = "Model") #+
#scale_x_continuous(breaks = seq(from = 1, to = 30, by = 2))

pdf(paste0(IncludeFigHere,"/PlotGenes.pdf"), 6,5 )
print(PlotGenes)
dev.off()

pdf(paste0(IncludeFigHere,"/PlotThree.pdf"), 15,5 )
print(PlotClinical + PlotImaging + PlotGenes)
dev.off()

Perf2 <- Perf %>%
  ggplot(aes(x=Workflow_Rank, y = mean, shape = Model_Type, color = Recipe)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-std_err, ymax = mean+std_err)) +
  theme_minimal()+
  #scale_colour_manual(values = cols) +
  labs(title = "Performance Comparison of Models", x = "Model Rank", y = "RMSE", color = "Data", shape = "Model") #+
#scale_x_continuous(breaks = seq(from = 1, to = 30, by = 2))

pdf(paste0(IncludeFigHere,"/Performance3.pdf"),7,6)
print(Perf2)
dev.off()

#stopCluster()

#################

vi_rnks2$var <-as.factor(vi_rnks2$var)

MostImpVars <- vi_rnks2 %>%
  filter(rnk < Threshold) %>%
  group_by(var) %>%
  count()


#######################

#Properly report on performance metrics - AIC and BIC 

Perf <- collect_metrics(race_results) 

save(Perf, file = paste0(IncludeFigHere,"/Prediction.RData"))

############################


#setwd("/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/Magda/Oct_20231028_2032") #this was enabled?
 
#load and get extra values - plus mediation analysis 

#install.packages("mediation")

#library(mediation)

#load("Prediction.RData")
#load("Modeling2data3mods.RData")
#load("ModelingWorkflows.RData")

#############################

Perf <- collect_predictions(race_results) 

#AIC and BIC- calculated for parametric models not RF , XGBoost etc. 


# can select different grids (JAMES)

#my_models <- 
#  workflow_set(
#    preproc = list(ibm_rec),
#    models = list(glmnet = log_spec,  xgbTree = xgb_spec),
#    cross = TRUE
#  ) %>%
#  # add custom grid 
#  option_add(grid = xgbTree_grid, id = "recipe_xgbTree") %>%
#  option_add(grid = glmnet_grid, id = "recipe_glmnet") 
#
#my_models


# optimize based on this? - 

#f1 <- function(data, lev = NULL, model = NULL) {
#  f1_val <- MLmetrics::F1_Score(y_pred = data$pred,
#                                y_true = data$obs,
#                                positive = lev[1])
#  c(F1 = f1_val)
#}

#ibm_metrics <- metric_set(bal_accuracy, roc_auc, yardstick::sensitivity, yardstick::specificity, yardstick::precision, f_meas)

#https://www.youtube.com/watch?v=2OfTEakSFXQ
###########

library(stacks)
#
ens <- stacks() %>%
  add_candidates(race_results) %>%
  blend_predictions() %>%
  fit_members()
#
##collect_parameters(ens,race_results )
#
#aa <- data_tst %>% 
#  rename( "Label" = "GeneralCog")
#
### why not all?
#Predict <- predict(ens, data_tst, outcomes = TRUE, members = TRUE) %>%
#  add_column(GeneralCog = data_tst$GeneralCog)
#
#PerformanceTest <- map_dfr(Predict, mae, truth = GeneralCog, data = Predict) %>%
#  mutate(member = colnames(Predict))


#
#member_preds <- 
#  tree_frogs_test %>%
#  select(latency) %>%
#  bind_cols(predict(tree_frogs_model_st, tree_frogs_test, members = TRUE))

save(ens,   file = paste0(IncludeFigHere,"/stacks.RData"))

pdf(paste0(IncludeFigHere,"/Stacks1.pdf"), 4,5)
autoplot(ens) + theme_bw()
dev.off()

pdf(paste0(IncludeFigHere,"/Stacks2.pdf") ,4,5)
autoplot(ens, type = "weights") + theme_bw()
dev.off()



#library(tidyposterior)
#
#rmse_set <- race_results() %>%
#  perf_mod(
#    iter = 5000, 
#    chains = 10. 
#    cores = 10, 
#    seed = 2, 
#    refresh = 0
#  )
#
#rmse_res %>%
#  autoplot(type = "ROPE", size = 0.15)

mae_stats <- collect_metrics(race_results, summarize = FALSE) %>%
  filter(.metric == "mae") %>%
  separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") 


PerfMAE <- list()

for (i in unique(mae_stats$Model_Type)){
  
  full_model <- filter(mae_stats, 
                       Recipe == "AllRec2J", Model_Type == i) %>%
    dplyr::select(id, full_model= .estimate)
  
  rest <- filter(mae_stats, 
                 Recipe != "AllRec2J", Model_Type == i) %>%
    dplyr::select(id, .estimate, Recipe)
  
  
  PerfMAE[[i]] <- left_join(rest, full_model, by = 'id') %>%
    mutate(loss = .estimate - full_model) %>%
    mutate(term =  Recipe) %>%
    nest_by(term) %>%
    mutate(stats = list(t.test(data$loss ~ 1, data = data))) %>%
    summarise(tidy(stats), .groups = "drop") %>%
    ggplot(aes(x= estimate, y = reorder(term , estimate))) +
    geom_point() + 
    geom_vline(xintercept = 0, lty = 2, col = "red") + 
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = .25) + 
    labs(x = "loss of MAE", y = NULL) + theme_bw() + ggtitle(paste0("Recipe ", i))
  
  
}

PerfMAE[["xgboost"]] <- NULL

pdf(paste0(IncludeFigHere, "/PerfromanceMAE.pdf"), 10, 7)
print(wrap_plots(PerfMAE))
dev.off()

# would be good to have tidyposterior to measure which one best in general - this one best in particular 


############### mediaton 

SelectFeatures <- result %>%
  bind_rows() %>%
  filter(rnk < 3) 

SelectFeatures <- unique(SelectFeatures$var) %>% 
  setdiff(c("Clinical_EducationMeasurement", "Clinical_Sex_X0", "e2e3_X0","rs17046359_G_X2"))



#https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171


ToMediate <- AllRec2J %>%
  dplyr::select(c(SelectFeatures), GeneralCog, Clinical_EducationMeasurement) #%>%
# rename(dv = GeneralCog,
#       iv = Clinical_EducationMeasurement)



#install.packages("psych")
library(psych)

# has a cool comorbidity function? 
#HAVE TO ZERO CENTRE: glbwarmc <-data.frame(scale(globalWarm,scale=FALSE))

#registerDoSEQ() 

mod.glb <- list()
Params <- list()

SelectFeatures <- setdiff(SelectFeatures, c("rs118186402_A_X2", "rs77433227_T_X2", "rs141892450_A_X2", "rs73353563_T_X2")) #inside calculation - Error in dimnames(x) <- dn : 
#length of 'dimnames' [2] not equal to array extent
for (i in SelectFeatures) {
  
  print(i)
  
  
  ToMediate2 <- ToMediate %>% 
    dplyr::select(GeneralCog, Clinical_EducationMeasurement, i)
  
  #debugonce(psych::mediate)
  mod.glb[[i]] <-  mediate(y = "GeneralCog", x = "Clinical_EducationMeasurement", m= eval(i), data = ToMediate2, plot = TRUE)
  
  Params[[i]] <- data.frame(Mean = mod.glb[[i]][["boot"]]$mean, Sd = mod.glb[[i]][["boot"]]$sd, ci.low = mod.glb[[i]][["boot"]]$ci[1,], 
                            ci.high = mod.glb[[i]][["boot"]]$ci[2,])
  
  
}


print("Mediation Figures")

pdf(paste0(IncludeFigHere, "/MediateDiagramAll.pdf"), 10, 5)

for (i in SelectFeatures){
  
  print(mediate.diagram(mod.glb[[i]]))
  
}

dev.off()

print("Mediation Figures Parameters")

ParamsPlot <-   Params %>%
  bind_rows(.,.id = 'Features') %>%
  mutate(Type = case_when(startsWith(Features, "WM") | startsWith(Features, "c") ~ "Networks", 
                          startsWith(Features, "rs") | startsWith(Features, "e") ~ "Genes", 
                          Features %in% c("ClinicalAge","Clinical_Sex_X1")  ~ "Demographics", 
                          Features == "Clinical_EducationMeasurement" ~ "Education",
                          TRUE ~ "Cognition")) %>%
  ggplot(aes(x=Features, y = `Mean.1`,color = Type)) + #Before Mean
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_point() +
  geom_errorbar(aes(ymin = `Mean.1` - `Sd.1`, ymax = `Mean.1` + `Sd.1`)) + theme_bw() + coord_flip()  

pdf(paste0(IncludeFigHere, "/MediationPerf.pdf"), 10, 5)
print(ParamsPlot)
dev.off()

#####################

#Coeffs <- list()

#library(mediation)
#
#
#fit.totaleffect <- lm(GeneralCog~Clinical_EducationMeasurement,ToMediate) 

#Coeffs[[2]] <- data.frame(Coeffs = fit.totaleffect$coefficients[["Clinical_EducationMeasurement"]])  %>% 
#  add_column(Type = "Second")

############
##########################

############################# Group features together? PCA 

#library(PCAtools) 

