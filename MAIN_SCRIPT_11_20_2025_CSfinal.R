## THE MAIN PURPOSE OF THIS SCRIPT IS TO PUT EVERYTHING INTO ONE CONSOLIDATED
## PLACE TO RUN THE CODE FOR THE ANALYSIS
#This is the code for:
#"Predictive Models of Seed Yield in Intermediate Wheatgrass (Thinopyrum intermedium): The Added Value of a Time Series of Multispectral Imagery in Semiarid Environments". 
# Claire Spickermann, Efrain Duarte, Steve Larson, Zayne Maughan and Alexander Hernandez





######### Libraries Needed
library(dplyr)
library(caret)
library(Metrics)
library(valmetrics)
library(creditmodel)
library(e1071)
library(doParallel)
library(tidyplots)
library(RColorBrewer)
library(readr)
library(stringr)
library(tidyr)
library(ggpubr)
library(randomForest)
library(reshape2)
########################### Internal FUNCTIONS ##################

##### adjust_cols is a function that transforms data into a better format
adjust_cols <- function(x) {
  x <- x %>%
    mutate(
      year = str_extract(colnames(.)[7], "[[:digit:]]{4}"),
      month = str_extract(
        str_extract(colnames(.)[7], "_[[:digit:]]{2}_"),
        "[[:digit:]]{2}"
      ),
      date = str_extract(colnames(.)[7], "\\d{4}_\\d{2}_\\d{2}")
    ) %>%
    rename_with(~ ifelse(str_detect(.x, "BGI2[[:digit:]]{4}"),
                         str_extract(., "[[:alpha:]]+[[:digit:]]"),
                         str_extract(.x, "[[:alpha:]]+")
    ))

  x
}


# metrics is a function that returns a named vector of metrics
### to evaluate models

metrics <- function(preds, actual, data) {
  rss <- sum((actual - preds)^2)
  tss <- sum((actual - mean(actual))^2)
  r_squared <- 1 - (rss / tss)
  n <- nrow(data)
  p <- ncol(data) - 1

  metric_list <- c(
    "R_squared" = r_squared,
    "adj_R_Squared" = 1 - ((1 - r_squared) * (n - 1)) / (n - p - 1),
    "mse" = mse(actual, preds),
    "rmse" = rmse(actual, preds),
    "mae" = mae(actual, preds),
    "lccc" = valmetrics::lccc(actual, preds)
  )
  metric_list
}

# this function is provided to eliminate large code sections by allowing a function
##  to be called
graph_type <- function(model, test_data, title = "Title",
                       graph_type = "VarImp", scale_adj = 0.1){
  #Get the predictions
  preds_model <- predict(model, test_data)

  #put the predictions and the response in one data frame for ease of graphics
  p_o_graph <- bind_cols(preds_model, test_data$adj_response) |>
    rename(
      preds = "...1",
      actual = "...2"
    ) |>
    # include the residuals for other graphics
    mutate(
      residuals = actual - preds
    )

  #needed for the preds v obs graph
  model_metrics <- metrics(preds_model, test_data$adj_response, test_data)

  #if (graph_type == "VarImp"){

    # variable importance graphs (from Zayne's master script)
    #var_imp_data <- varImp(model)$importance

    #var_imp_data <- data.frame(
      #variables = row.names(var_imp_data),
      #importance = var_imp_data$Overall
    #) |> dplyr::arrange(importance)

    #var_imp_data$variables <- factor(var_imp_data$variables,
                                     #levels = var_imp_data$variables)

    #variable importance graphs (From Efrain's august 15 script)
  if (graph_type == "VarImp") {
      
      var_imp_data <- varImp(model)$importance |> 
        tibble::rownames_to_column("variables") |> 
        arrange(Overall) |>
        rename(importance = Overall) #|> # remove hashtag to get the reduced variable importance chart
        #slice_tail(n = 10)
      
      var_imp_data$variables <- factor(var_imp_data$variables,
                                       levels = var_imp_data$variables)
    graph <-  var_imp_data |>
      tidyplot(x = importance, y = variables, color = variables) |>
      add_mean_bar() |>
      adjust_colors(new_colors = brewer.pal(7, "Blues")) |>
      adjust_size(NA, NA) |>
      adjust_font(fontsize = 12) |>
      adjust_x_axis_title(title = "Mean Decrease in Accuracy", fontsize = 12) |>
      # adjust_x_axis(limits = c(0, 20000000))|>
      adjust_y_axis_title(title = "Vegetation Index", fontsize = 12) |>
      add_title(title) |>
      remove_legend() |>
      add(
        ggplot2::theme(
          plot.title   = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.x  = element_text(size = 10),
          axis.text.y  = element_text(size = 10)
        )
      )



  } else if (graph_type == "p_v_o"){
    ###coordenadas para el modelo de regresión anual unmask
    xpos <- max(p_o_graph$actual) - 0.05 * sd(p_o_graph$actual)
    ypos_r2 <- min(p_o_graph$preds) + 0.45 * sd(p_o_graph$preds)
    ypos_ci <- min(p_o_graph$preds) + 0.10 * sd(p_o_graph$preds)

    graph <- p_o_graph |>
      tidyplot(x = actual, y = preds) |>
      add_data_points(color = "#9ecae1") |>
      add_curve_fit(se = TRUE, color = "#3182bd", fill = "#9ecae1") |>
      add_annotation_text(
        str_c("R²: ", formatC(model_metrics["R_squared"], digits = 3, format = "f")),
        x = xpos,
        y = ypos_r2,
        fontsize = 12
      ) |>
      add_annotation_text(
        "CI: 95%",
        x = xpos,
        y = ypos_ci,
        fontsize = 11
      ) |>
      adjust_size(NA, NA) |>
      adjust_font(fontsize = 12) |>
      adjust_x_axis_title("Observed", fontsize = 12) |>
      adjust_y_axis_title("Predicted", fontsize = 12) |>
      add_title(title) |>
      add(
        ggplot2::theme(
          plot.title   = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.x  = element_text(size = 10),
          axis.text.y  = element_text(size = 10)
        )
      )
  } else if (graph_type == "residuals") {
    
    graph <- p_o_graph |>
      tidyplot(x = actual, y = residuals) |>
      add_data_points(color = "#9ecae1") |>
      add_curve_fit(se = FALSE, color = "#3182bd") |>
      adjust_size(NA, NA) |>
      adjust_font(fontsize = 12) |>
      adjust_x_axis_title("Observed", fontsize = 12) |>
      adjust_y_axis_title("Residuals (Observed - Predicted)", fontsize = 12) |>
      add_title(title) |>
      add(
        ggplot2::theme(
          plot.title   = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.x  = element_text(size = 10),
          axis.text.y  = element_text(size = 10)
        )
      )
  }
  
  return(graph)
}
######################## Read in the data #############################
# READ IN the MASKED DATA

data_path <- "data/masked_data"
files <- list.files(data_path, pattern = "*.csv", full.names = TRUE)
data <- lapply(files, read_csv)
data <- lapply(data, adjust_cols) |> bind_rows()


## READ IN THE DATA FOR THE UNMASKED DATA
## create a path of the files for future use
unmasked_data_path <- "data/unmasked_data"
unmasked_files <- list.files(unmasked_data_path,
  pattern = "*.csv",
  full.names = TRUE
)
# read in each file
unmasked_data <- lapply(unmasked_files, read_csv)
# convert the data itself and then combine all items in the list
unmasked_data <- lapply(unmasked_data, adjust_cols) |> bind_rows()

unmasked_data <- unmasked_data[,c(1:10,12:38)] #Remove SCI
unmasked_data <- unmasked_data[,c(1:30,32:37)] #remove BGI2
unmasked_data <- unmasked_data[,c(1:28,30:36)] #remove GR
unmasked_data <- unmasked_data[,c(1:33,35)]#remove GLA
############################### AREA FOR MASKED PLOTS ########################

## GET THE AREA FOR THE MASKED PLOTS
area_data_path <- "data/masked_area"
area_dat_files <- list.files(area_data_path,
                             pattern = "*csv",
                             full.names = TRUE)


# READ IN THE DATA AND EXTRACT THE DATE TO COMBINE THIS WITH THE OTHER DATA
area_dat <- lapply(area_dat_files, function(file) {
  read_csv(file) |>
    mutate(
      year = str_extract(str_extract(file, "20\\d{2}_"), "\\d{4}"),
      month = str_extract(file, "(?<=_)\\d{2}"),
      day = str_extract(file, "(?<=_)\\d{4}") |> str_sub(3, 4),
      date = str_c(year, month, day, sep = "_")
    )
}) |>
  bind_rows() |>
  mutate(area = ifelse(!is.na(area),
    area,
    ifelse(is.na(area_calculation),
      `area(m^2)`,
      area_calculation
    )
  )) |>
  select(plant_id, area, id, year, month, date)

# NEED THIS TO ENSURE COHESIVE JOININGS
Plant_id <- read_csv("data/Plant_ID_key.csv")

########################## GET THE RESPONSE or YIELD ######################
#### EACH YEAR HAD A DIFFERENT YIELD THAT WAS RECORDED. THIS CHUNK READS IN
####   AND EXTRACTS THE NEEDED INFORMATION
###### TO COMBINE IT WITH OTHER DATA LATER


#### FOCUS ON THE RESPONSE VARIABLE
response_2024 <- read_csv("data/gsd_plot2024_Richmond_editted.csv") |>
  mutate(
    plant_id = plot,
    year = "2024",
    gsdplot = replace_na(gsdplot, "0")
  ) |>
  mutate(
    gsdplot = ifelse(gsdplot == "na", "0", gsdplot)
  ) |>
  select(plant_id, gsdplot, year)
############ 2023
response_2023 <- read_csv("data/Richmond Weights 2023-27OCT2023.csv") |>
  rename(weight = `Weight (g)`) |>
  slice_tail(n = -1) |>
  mutate(
    plant_id = ID,
    year = "2023"
  ) |>
  mutate(
    plant_id = ifelse(plant_id == "2RITI21649", "20RITI21649", plant_id)
  ) |>
  select(weight, plant_id, year)

############## 2022
response_2022 <- read_csv("data/Kernza_2022_Richmond_Seed_Weights.csv") |>
  mutate(
    plant_id = entity_id,
    year = "2022"
  ) |>
  select(phenotype_value, plant_id, year)

###### COMBINE THE YIELDS INTO ONE SINGULAR TABLE SO IT CAN BE COMBINED LATER
response <- response_2024 |>
  inner_join(response_2023, by = c("plant_id")) |>
  inner_join(response_2022, by = c("plant_id")) |>
  pivot_wider(values_from = phenotype_value, names_from = year) |>
  pivot_wider(values_from = gsdplot, names_from = year.x) |>
  pivot_wider(values_from = weight, names_from = year.y) |>
  rename(
    response_2024 = `2024`,
    response_2023 = `2023`,
    response_2022 = `2022`
  )

response <-  response |>
  inner_join(Plant_id, by = c("plant_id")) |>
  select(
    plant_id, response_2022,
    response_2023, response_2024, id
  )


###############################################################################

# READ IN THE PLANT ID AND PLOTS THAT ARE NEEDED TO FILTER OUT OF THE DATA
alkar <- read_csv("data/Alkar_locations.csv") |> drop_na()
Plant_id <- read_csv("data/Plant_ID_key.csv")
plant_list <- read_csv("data/2024_plot_trimmed_list.csv")

# ALKAR DID NOT HAVE AN ID COLUMN SO THERE WAS NO WAY TO REFERENCE it with the
# DATA

#alkar_id <- alkar |>
  #dplyr::inner_join(plant_id, by = c("plant_id")) |>
  #select(id)

#Editting code to have it combine by the row range not id
alkar_id <- alkar |>
  dplyr::inner_join(Plant_id, by = c("plant_id")) 
  
plant_list <- plant_list |>
  dplyr::inner_join(Plant_id, by = c("Row_Range")) 
################################################################################
combined_data <- data |>
  inner_join(response, by = c("id" = "id")) |>
  mutate(
    response = as.numeric(ifelse(year == "2024", response_2024,
      ifelse(year == "2023", response_2023, response_2022)
    ))
  ) |>
  select(-c(response_2024, response_2023, response_2022))


masked_dataset <- combined_data |>
  inner_join(area_dat, by = c("plant_id", "date")) |>
  ## CONVERT TO KG/HA
  mutate(adj_response = (response * 10 / 3)) |>
  dplyr::filter(!(plant_id %in% unique(plant_list$plant_id)) &
                  !(plant_id %in% alkar_id$plant_id)) |>
  drop_na()
##### UNMASKED SET
unmasked_dataset <- unmasked_data |>
  mutate(
    area = (right - left) * (top - bottom)
  ) |>
  inner_join(response, by = c("id")) |>
  mutate(
    response = as.numeric(ifelse(year == "2024",
      response_2024,
      ifelse(year == "2023",
        response_2023, response_2022
      )
    ))
  ) |>
  select(-c(response_2024, response_2023, response_2022)) |>
  ### CONVERT TO KG/HA
  mutate(adj_response = (response * 10 / 3)) |>
  dplyr::filter(!(plant_id %in% unique(plant_list$plant_id)) &
                  !(plant_id %in% alkar_id$plant_id)) |>
  drop_na()


# adjust the plots to get it to match up


# FILTER OUT THE ROWS OF THE UNDESIRED PLANTS whether from Control group or were
# harvested at an earlier time than the other plants
#dataset <- data |>
  #dplyr::filter(!(id %in% unique(plant_list$id)) &
                  #!(id %in% alkar_id))

############################################################################
#SEED YEILD GRAPHS AND FIGURES
tidyplot(unmasked_dataset, x = year, y = adj_response ) |>
  add_boxplot() |> adjust_size(NA, NA) |>
  adjust_title(title = "Boxplot of Yield per year (kg/ha) Unmasked",
               fontsize = 10) |>
  adjust_y_axis_title(title = "Yield (kg/ha)", fontsize = 8) |>
  adjust_x_axis_title("Year", fontsize = 8)

summary(unmasked_dataset$adj_response)

tidyplot(masked_dataset, x = year.x, y = adj_response ) |>
  add_boxplot() |> adjust_size(NA, NA) |>
  adjust_title(title = "Boxplot of Yield per year (kg/ha) (Masked)",
               fontsize = 10) |>
  adjust_y_axis_title(title = "Yield (kg/ha)", fontsize = 8) |>
  adjust_x_axis_title("Year", fontsize = 8)

summary(masked_dataset$adj_response)

############################################################################
############################# MODELING #####################################

################################## Pre-Harvest Flights All indices
unmasked_dataset |>
  dplyr::filter( date %in% c("2022_08_07", "2023_08_07", "2024_07_10" )) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))  |>
  train_test_split(prop = 0.7)  -> unmasks_harvest

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)

# FIT THE MODEL

# TO SEE THE MODEL BEFORE THESE SELECT PREDICTORS: COMMMENT OUT THE |> AND THE
#  select() statement
rf_model_harvest <- train(adj_response ~ .,
                               data = unmasks_harvest$train, 
                               method = "rf",
                               trControl = repeat_oob,
                               metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_rf_harvest_globalpre <- predict(rf_model_harvest, newdata = unmasks_harvest$test)
unmask_harvest <- metrics(preds_rf_harvest_globalpre, unmasks_harvest$test$adj_response,
                          unmasks_harvest$test)

residuals_harv <- preds_rf_harvest_globalpre - unmasks_harvest$test$adj_response

summ_residuals_harv <- data.frame(mean = mean(residuals_harv), sd = sd(residuals_harv),
                                       max = max(residuals_harv), min = min(residuals_harv),
                                       range = max(residuals_harv) - min(residuals_harv))
############################### variable importance
p1 <- graph_type(rf_model_harvest, unmasks_harvest$test,
           title = "Pre-Harvest Global Model: Variable Importance", graph_type = "VarImp")
print(p1)
################################# preds vs obs
p2 <- graph_type(rf_model_harvest, unmasks_harvest$test,
           title = "Pre-Harvest Global Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
p3 <- graph_type(rf_model_harvest, unmasks_harvest$test,
           title = "Pre-Harvest Global Model",
           graph_type = "residuals")


ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/Pre_Harvest_all_VariableImportance_unmask_wo_SCI_fixed_alkar_reduced.jpg", plot = p1, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/Pre_Harvest_all_Observed_vs_Predicted_unmask_wo_SCI_fixed_alkar.jpg", plot = p2, device = "jpeg", dpi = 300,
      # width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/Pre_Harvest_all_Residuals_unmask_wo_SCI_fixed_alkar.jpg", plot = p3, device = "jpeg", dpi = 300,
      # width = 180, height = 120, units = "mm")

rsquared_table_preharvest <- bind_rows(unmask_harvest)


residuals_table_preharvest<- bind_rows(summ_residuals_harv)


residuals_table_preharvest
write.csv(rsquared_table_preharvest, file = "./residuals/all_indices_wo_SCI_BGI2_GR_GLA/rsquared_tables_preharvest_all_fixed_alkar.csv")
write.csv(residuals_table_preharvest, file = "./residuals/all_indices_wo_SCI_BGI2_GR_GLA/residuals_tables_preharvest_all_fixed_alkar.csv")


################################## July Flights
unmasked_dataset |>
  dplyr::filter( date %in% c("2022_07_11", "2023_07_05", "2024_07_10" )) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))  |>
  train_test_split(prop = 0.7)  -> unmasks_harvest

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)

# FIT THE MODEL

# TO SEE THE MODEL BEFORE THESE SELECT PREDICTORS: COMMMENT OUT THE |> AND THE
#  select() statement
rf_model_harvest_July <- train(adj_response ~ .,
                          data = unmasks_harvest$train, # |>
                            #select(GNDVI, CVI, CIG,NDVI, DVI, RVI, TVI, IHS, BI,
                                   #adj_response),
                          method = "rf",
                          trControl = repeat_oob,
                          metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_rf_harvest_globalJuly <- predict(rf_model_harvest_July, newdata = unmasks_harvest$test)
unmask_harvest_July <- metrics(preds_rf_harvest_globalJuly, unmasks_harvest$test$adj_response,
                          unmasks_harvest$test)

residuals_harv <- preds_rf_harvest_globalJuly - unmasks_harvest$test$adj_response

summ_residuals_harv_July <- data.frame(mean = mean(residuals_harv), sd = sd(residuals_harv),
                          max = max(residuals_harv), min = min(residuals_harv),
                           range = max(residuals_harv) - min(residuals_harv))
############################### variable importance
p7 <- graph_type(rf_model_harvest_July, unmasks_harvest$test,
           title = "July Global Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
p8 <- graph_type(rf_model_harvest_July, unmasks_harvest$test,
           title = "July Global Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals 
p9 <- graph_type(rf_model_harvest_July, unmasks_harvest$test,
           title = "July Global Model",
           graph_type = "residuals")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July_all_VariableImportance_unmask_wo_SCI_fixed_alkar_reduced.jpg", plot = p7, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July_all_Observed_vs_Predicted_unmask_wo_SCI_fixed_alkar.jpg", plot = p8, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July_all_Residuals_unmask_wo_SCI_fixed_alkar.jpg", plot = p9, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

################################## June Flights
unmasked_dataset |>
  dplyr::filter( date %in% c("2022_06_21", "2023_06_16", "2024_06_10" )) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))  |>
  train_test_split(prop = 0.7)  -> unmasks_harvest


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)

# FIT THE MODEL

# TO SEE THE MODEL BEFORE THESE SELECT PREDICTORS: COMMMENT OUT THE |> AND THE
#  select() statement
rf_model_harvest_June <- train(adj_response ~ .,
                          data = unmasks_harvest$train, 
                          method = "rf",
                          trControl = repeat_oob,
                          metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_rf_harvest_globalJune <- predict(rf_model_harvest_June, newdata = unmasks_harvest$test)
unmask_harvest_June <- metrics(preds_rf_harvest_globalJune, unmasks_harvest$test$adj_response,
                          unmasks_harvest$test)

residuals_harv <- preds_rf_harvest_globalJune - unmasks_harvest$test$adj_response

summ_residuals_harv_June <- data.frame(mean = mean(residuals_harv), sd = sd(residuals_harv),
                                  max = max(residuals_harv), min = min(residuals_harv),
                                  range = max(residuals_harv) - min(residuals_harv))
############################### variable importance
p13 <- graph_type(rf_model_harvest_June, unmasks_harvest$test,
           title = "June Global Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
p14 <- graph_type(rf_model_harvest_June, unmasks_harvest$test,
           title = "June Global Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
p15 <- graph_type(rf_model_harvest_June, unmasks_harvest$test,
           title = "June Global Model",
           graph_type = "residuals")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June_all_VariableImportance_unmask_wo_SCI_fixed_alkar_reduced.jpg", plot = p13, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June_all_Observed_vs_Predicted_unmask_wo_SCI_fixed_alkar.jpg", plot = p14, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June_all_Residuals_unmask_wo_SCI_fixed_alkar.jpg", plot = p15, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")


################################## May Flights
unmasked_dataset |>
  dplyr::filter( date %in% c("2022_05_13", "2023_05_19", "2024_05_20" )) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))  |>
  train_test_split(prop = 0.7)  -> unmasks_harvest

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)

# FIT THE MODEL

# TO SEE THE MODEL BEFORE THESE SELECT PREDICTORS: COMMMENT OUT THE |> AND THE
#  select() statement
rf_model_harvest_May <- train(adj_response ~ .,
                               data = unmasks_harvest$train,
                               method = "rf",
                               trControl = repeat_oob,
                               metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_rf_harvest_globalMay <- predict(rf_model_harvest_May, newdata = unmasks_harvest$test)
unmask_harvest_May <- metrics(preds_rf_harvest_globalMay, unmasks_harvest$test$adj_response,
                          unmasks_harvest$test)

residuals_harv <- preds_rf_harvest_globalMay - unmasks_harvest$test$adj_response

summ_residuals_harv_May <- data.frame(mean = mean(residuals_harv), sd = sd(residuals_harv),
                                       max = max(residuals_harv), min = min(residuals_harv),
                                       range = max(residuals_harv) - min(residuals_harv))
############################### variable importance
p19 <- graph_type(rf_model_harvest_May, unmasks_harvest$test,
           title = "May Global Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
p20 <- graph_type(rf_model_harvest_May, unmasks_harvest$test,
           title = "May Global Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
p21 <- graph_type(rf_model_harvest_May, unmasks_harvest$test,
           title = "May Global Model",
           graph_type = "residuals")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May_all_VariableImportance_unmask_wo_SCI_fixed_alkar_reduced.jpg", plot = p19, device = "jpeg", dpi = 300, 
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May_all_Observed_vs_Predicted_unmask_wo_fixed_alkar.jpg", plot = p20, device = "jpeg", dpi = 300, 
      width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May_all_Residuals_unmask_wo_SCI_fixed_alkar.jpg", plot = p21, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")


metrics_months <- bind_rows(unmask_harvest_May, unmask_harvest_June,  unmask_harvest_July) |>
  bind_cols(c("Unmask May",  "Unmask June", "Unmask July")) |> rename("Type" = "...7") |>
  relocate(Type)

residuals_table_JuneJuly <- bind_rows(summ_residuals_harv_May, summ_residuals_harv_June, summ_residuals_harv_July) |>
  bind_cols(c("Unmasked May", "Unmasked June", 
               "Unmasked July")) |> rename("Model" = "...6") |>
  relocate(Model)


residuals_table_JuneJuly
write.csv(metrics_months, file = "./residuals/all_indices_wo_SCI_BGI2_GR_GLA/rsquared_tables_months_all_fixed_alkar.csv")
write.csv(residuals_table_JuneJuly, file = "./residuals/all_indices_wo_SCI_BGI2_GR_GLA/residuals_tables_months_all_fixed_alkar.csv")

################################################################################
#########Recreating the individual R^2 values graph Alex had ###################

#########################May 2022###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2022_05_13")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv22


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_May22 <- train(adj_response ~ .,
                    data = test_indiv22$train,
                    method = "rf",
                    trControl = repeat_oob,
                    metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_May22 <- predict(rf_harv_May22, newdata = test_indiv22$test)
rf_metrics_May22 <- metrics(preds_test_May22, test_indiv22$test$adj_response,
                      test_indiv22$test)
residuals_harv <- preds_test_May22 - test_indiv22$test$adj_response

summ_residuals_harv_May22 <- data.frame(mean = mean(residuals_harv),
                                    sd = sd(residuals_harv),
                                    max = max(residuals_harv),  min = min(residuals_harv),
                                    range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s1 <- graph_type(rf_harv_May22, test_indiv22$test,
           title = "May 2022 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s2 <- graph_type(rf_harv_May22, test_indiv22$test,
           title = "May 2022 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s3 <- graph_type(rf_harv_May22, test_indiv22$test, title = "May 2022 Model",
           graph_type = "residuals")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May22_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar_reduced.jpg", plot = s1, device = "jpeg", dpi = 300,
         width = 180, height = 120, units = "mm")
       
ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May22_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s2, device = "jpeg", dpi = 300,
        width = 180, height = 120, units = "mm")
       
ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May22_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s3, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")
#########################June 2022###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2022_06_21")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv22


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_Jun22 <- train(adj_response ~ .,
                       data = test_indiv22$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_Jun22 <- predict(rf_harv_Jun22, newdata = test_indiv22$test)
rf_metrics_Jun22 <- metrics(preds_test_Jun22, test_indiv22$test$adj_response,
                            test_indiv22$test)
residuals_harv <- preds_test_Jun22 - test_indiv22$test$adj_response

summ_residuals_harv_Jun22 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s4 <- graph_type(rf_harv_Jun22, test_indiv22$test,
           title = "June 2022 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s5 <- graph_type(rf_harv_Jun22, test_indiv22$test,
           title = "June 2022 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s6 <- graph_type(rf_harv_Jun22, test_indiv22$test, title = "June 2022 Model",
           graph_type = "residuals")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June22_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s4, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June22_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s5, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June22_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s6, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")
#########################July 2022###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2022_07_11")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv22


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_July22 <- train(adj_response ~ .,
                       data = test_indiv22$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_July22 <- predict(rf_harv_July22, newdata = test_indiv22$test)
rf_metrics_July22 <- metrics(preds_test_July22, test_indiv22$test$adj_response,
                            test_indiv22$test)
residuals_harv <- preds_test_July22 - test_indiv22$test$adj_response

summ_residuals_harv_July22 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s7 <- graph_type(rf_harv_July22, test_indiv22$test,
           title = "July 2022 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s8 <- graph_type(rf_harv_July22, test_indiv22$test,
           title = "July 2022 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s9 <- graph_type(rf_harv_July22, test_indiv22$test, title = "July 2022 Model",
           graph_type = "residuals")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July22_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s7, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July22_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s8, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July22_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s9, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")
#########################August 2022###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2022_08_07")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv22


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_Aug22 <- train(adj_response ~ .,
                       data = test_indiv22$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_Aug22 <- predict(rf_harv_Aug22, newdata = test_indiv22$test)
rf_metrics_Aug22 <- metrics(preds_test_Aug22, test_indiv22$test$adj_response,
                            test_indiv22$test)
residuals_harv <- preds_test_Aug22 - test_indiv22$test$adj_response

summ_residuals_harv_Aug22 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s10 <- graph_type(rf_harv_Aug22, test_indiv22$test,
           title = "August 2022 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s11 <- graph_type(rf_harv_Aug22, test_indiv22$test,
           title = "August 2022 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s12 <- graph_type(rf_harv_Aug22, test_indiv22$test, title = "August 2022 Model",
           graph_type = "residuals")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/August22_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s10, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/August22_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s11, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/August22_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s12, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")
#########################April 2023###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2023_05_03")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv23


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_Apr23 <- train(adj_response ~ .,
                    data = test_indiv23$train,
                    method = "rf",
                    trControl = repeat_oob,
                    metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_Apr23 <- predict(rf_harv_Apr23, newdata = test_indiv23$test)
rf_metrics_Apr23 <- metrics(preds_test_Apr23, test_indiv23$test$adj_response,
                      test_indiv23$test)
residuals_harv <- preds_test_Apr23 - test_indiv23$test$adj_response

summ_residuals_harv_Apr23 <- data.frame(mean = mean(residuals_harv),
                                    sd = sd(residuals_harv),
                                    max = max(residuals_harv),  min = min(residuals_harv),
                                    range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s13 <- graph_type(rf_harv_Apr23, test_indiv23$test,
           title = "April 2023 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s14 <- graph_type(rf_harv_Apr23, test_indiv23$test,
           title = "April 2023 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s15 <- graph_type(rf_harv_Apr23, test_indiv23$test, title = "April 2023 Model",
           graph_type = "residuals")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/April23_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s13, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/April23_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s14, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/April23_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s15, device = "jpeg", dpi = 300,
    #   width = 180, height = 120, units = "mm")
#########################May 2023###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2023_05_19")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv23


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_May23 <- train(adj_response ~ .,
                       data = test_indiv23$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_May23 <- predict(rf_harv_May23, newdata = test_indiv23$test)
rf_metrics_May23 <- metrics(preds_test_May23, test_indiv23$test$adj_response,
                            test_indiv23$test)
residuals_harv <- preds_test_May23 - test_indiv23$test$adj_response

summ_residuals_harv_May23 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s16 <- graph_type(rf_harv_May23, test_indiv23$test,
           title = "May 2023 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s17 <- graph_type(rf_harv_May23, test_indiv23$test,
           title = "May 2023 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s18 <- graph_type(rf_harv_May23, test_indiv23$test, title = "May 2023 Model",
           graph_type = "residuals")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May23_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s16, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May23_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s17, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May23_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s18, device = "jpeg", dpi = 300,
      # width = 180, height = 120, units = "mm")
#########################June 2023###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2023_06_16")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv23


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_June23 <- train(adj_response ~ .,
                       data = test_indiv23$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_June23 <- predict(rf_harv_June23, newdata = test_indiv23$test)
rf_metrics_June23 <- metrics(preds_test_June23, test_indiv23$test$adj_response,
                            test_indiv23$test)
residuals_harv <- preds_test_June23 - test_indiv23$test$adj_response

summ_residuals_harv_June23 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s19 <- graph_type(rf_harv_June23, test_indiv23$test,
           title = "June 2023 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s20 <- graph_type(rf_harv_June23, test_indiv23$test,
           title = "June 2023 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s21 <- graph_type(rf_harv_June23, test_indiv23$test, title = "June 2023 Model",
           graph_type = "residuals")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June23_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s19, device = "jpeg", dpi = 300,
    #   width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June23_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s20, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June23_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s21, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")
#########################July 2023###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2023_07_05")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv23


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_July23 <- train(adj_response ~ .,
                       data = test_indiv23$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_July23 <- predict(rf_harv_July23, newdata = test_indiv23$test)
rf_metrics_July23 <- metrics(preds_test_July23, test_indiv23$test$adj_response,
                            test_indiv23$test)
residuals_harv <- preds_test_July23 - test_indiv23$test$adj_response

summ_residuals_harv_July23 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s22 <- graph_type(rf_harv_July23, test_indiv23$test,
           title = "July 2023 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s23 <- graph_type(rf_harv_July23, test_indiv23$test,
           title = "July 2023 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s24 <- graph_type(rf_harv_July23, test_indiv23$test, title = "July 2023 Model",
           graph_type = "residuals")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July23_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s22, device = "jpeg", dpi = 300,
    #   width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July23_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s23, device = "jpeg", dpi = 300,
   #    width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July23_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s24, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")
#########################August 2023###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2023_08_07")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv23


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_Aug23 <- train(adj_response ~ .,
                       data = test_indiv23$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_Aug23 <- predict(rf_harv_Aug23, newdata = test_indiv23$test)
rf_metrics_Aug23 <- metrics(preds_test_Aug23, test_indiv23$test$adj_response,
                            test_indiv23$test)
residuals_harv <- preds_test_Aug23 - test_indiv23$test$adj_response

summ_residuals_harv_Aug23 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s25 <- graph_type(rf_harv_Aug23, test_indiv23$test,
           title = "August 2023 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s26 <- graph_type(rf_harv_Aug23, test_indiv23$test,
           title = "August 2023 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s27 <- graph_type(rf_harv_Aug23, test_indiv23$test, title = "August 2023 Model",
           graph_type = "residuals")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/August23_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s25, device = "jpeg", dpi = 300,
    #   width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/August23_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s26, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/August23_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s27, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")
#########################April 2024###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2024_04_24")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv24


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_Apr24 <- train(adj_response ~ .,
                    data = test_indiv24$train,
                    method = "rf",
                    trControl = repeat_oob,
                    metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_Apr24 <- predict(rf_harv_Apr24, newdata = test_indiv24$test)
rf_metrics_Apr24 <- metrics(preds_test_Apr24, test_indiv24$test$adj_response,
                      test_indiv24$test)
residuals_harv <- preds_test_Apr24 - test_indiv24$test$adj_response

summ_residuals_harv_Apr24 <- data.frame(mean = mean(residuals_harv),
                                    sd = sd(residuals_harv),
                                    max = max(residuals_harv),  min = min(residuals_harv),
                                    range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s28 <- graph_type(rf_harv_Apr24, test_indiv24$test,
           title = "April 2024 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s29 <- graph_type(rf_harv_Apr24, test_indiv24$test,
           title = "April 2024 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s30 <- graph_type(rf_harv_Apr24, test_indiv24$test, title = "April 2024 Model",
           graph_type = "residuals")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/April24_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s28, device = "jpeg", dpi = 300,
    #   width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/April24_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s29, device = "jpeg", dpi = 300,
    #   width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/April24_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s30, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")
#########################May 2024###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2024_05_20")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv24


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_May24 <- train(adj_response ~ .,
                       data = test_indiv24$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_May24 <- predict(rf_harv_May24, newdata = test_indiv24$test)
rf_metrics_May24 <- metrics(preds_test_May24, test_indiv24$test$adj_response,
                            test_indiv24$test)
residuals_harv <- preds_test_May24 - test_indiv24$test$adj_response

summ_residuals_harv_May24 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s31 <- graph_type(rf_harv_May24, test_indiv24$test,
           title = "May 2024 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s32 <- graph_type(rf_harv_May24, test_indiv24$test,
           title = "May 2024 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s33 <- graph_type(rf_harv_May24, test_indiv24$test, title = "May 2024 Model",
           graph_type = "residuals")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May24_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s31, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May24_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s32, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/May24_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s33, device = "jpeg", dpi = 300,
      # width = 180, height = 120, units = "mm")
#########################June 2024###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2024_06_10")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv24


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_June24 <- train(adj_response ~ .,
                       data = test_indiv24$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_June24 <- predict(rf_harv_June24, newdata = test_indiv24$test)
rf_metrics_June24 <- metrics(preds_test_June24, test_indiv24$test$adj_response,
                            test_indiv24$test)
residuals_harv <- preds_test_June24 - test_indiv24$test$adj_response

summ_residuals_harv_June24 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s34 <- graph_type(rf_harv_June24, test_indiv24$test,
           title = "June 2024 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s35 <- graph_type(rf_harv_June24, test_indiv24$test,
           title = "June 2024 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s36 <- graph_type(rf_harv_June24, test_indiv24$test, title = "June 2024 Model",
           graph_type = "residuals")

ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June24_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s34, device = "jpeg", dpi = 300,
       width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June24_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s35, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/June24_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s36, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")
#########################July 2024###############################################
unmasked_dataset |>
  dplyr::filter( date %in% c("2024_07_10")) |>
  select(-c(left, right, fid, top,
            bottom, response, plant_id,
            fid, id, year, month, date, area))|>
  train_test_split() -> test_indiv24


# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_oob <- trainControl(method = "oob", number = 10)


rf_harv_July24 <- train(adj_response ~ .,
                       data = test_indiv24$train,
                       method = "rf",
                       trControl = repeat_oob,
                       metric = "Rsquared"
)

# MAKE THE PREDICTIONS
preds_test_July24 <- predict(rf_harv_July24, newdata = test_indiv24$test)
rf_metrics_July24 <- metrics(preds_test_July24, test_indiv24$test$adj_response,
                            test_indiv24$test)
residuals_harv <- preds_test_July24 - test_indiv24$test$adj_response

summ_residuals_harv_July24 <- data.frame(mean = mean(residuals_harv),
                                        sd = sd(residuals_harv),
                                        max = max(residuals_harv),  min = min(residuals_harv),
                                        range = max(residuals_harv) - min(residuals_harv))
############################### variable importance

s37 <- graph_type(rf_harv_July24, test_indiv24$test,
           title = "July 2024 Model: Variable Importance", graph_type = "VarImp")

################################# preds vs obs
s38 <- graph_type(rf_harv_July24, test_indiv24$test,
           title = "July 2024 Model", graph_type = "p_v_o",
           scale_adj = 0.1)

#################################### Residuals
s39 <- graph_type(rf_harv_July24, test_indiv24$test, title = "July 2024 Model",
           graph_type = "residuals")

metrics_months <- bind_rows(rf_metrics_May22, rf_metrics_Jun22,  rf_metrics_July22, rf_metrics_Aug22, 
                            rf_metrics_Apr23, rf_metrics_May23, rf_metrics_June23, rf_metrics_July23, 
                            rf_metrics_Aug23, rf_metrics_Apr24, rf_metrics_May24, rf_metrics_June24, rf_metrics_July24) |>
  bind_cols(c("May 2022", "June 2022", "July 2022", "August 2022", "April 2023", "May 2023", 
              "June 2023", "July 2023", "August 2023", "April 2024", "May 2024", "June 2024",
              "July 2024")) |> rename("Type" = "...7") |>
  relocate(Type)

residuals_table_months <- bind_rows(summ_residuals_harv_May22, summ_residuals_harv_Jun22,
                             summ_residuals_harv_July22, summ_residuals_harv_Aug22,
                             summ_residuals_harv_Apr23, summ_residuals_harv_May23, 
                             summ_residuals_harv_June23, summ_residuals_harv_July23, 
                             summ_residuals_harv_Aug23,  summ_residuals_harv_Apr24,
                             summ_residuals_harv_May24, summ_residuals_harv_June24, 
                             summ_residuals_harv_July24) |>
  bind_cols(c("May 2022", "June 2022", "July 2022", "August 2022", "April 2023", "May 2023", 
              "June 2023", "July 2023", "August 2023", "April 2024", "May 2024", "June 2024",
              "July 2024")) |> rename("Model" = "...6") |>
  relocate(Model)


residuals_table_months

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July24_all_VariableImportance_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s37, device = "jpeg", dpi = 300,
      # width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July24_all_Observed_vs_Predicted_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s38, device = "jpeg", dpi = 300,
     #  width = 180, height = 120, units = "mm")

#ggsave("./plots/plots_fixed_alkar_no_SCI_BGI2_GR_GLA_fixed_axis/July24_all_Residuals_single_wo_SCI_notreduced_fixed_alkar.jpg", plot = s39, device = "jpeg", dpi = 300,
      # width = 180, height = 120, units = "mm")

#write.csv(metrics_months, file = "./residuals/all_indices_wo_SCI_BGI2_GR_GLA/rsquared_tables_months_all_single_flights_fixed_alkar.csv")
#write.csv(residuals_table_months, file = "./residuals/all_indices_wo_SCI_BGI2_GR_GLA/residuals_tables_months_all_single_flights_fixed_alkar.csv")





