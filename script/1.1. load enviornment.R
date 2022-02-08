# Corestem project - CA Analysis.proj
# 20220208 Sehwan Chun at Corestem, Inc.
# 1.1. load environment

#### 1. Library Loading ####
library("readxl")

#### 2. Files Loading ####

# Load Gene expression data 
expr_file = "./data/gene_count_matrix_211012.csv"
expr_file = read.csv(expr_file, stringsAsFactors = F, header = T)
expr_file = expr_file[complete.cases(expr_file),]

# Load Clinical data (SARA)
efficacy_file = "./data/CA_efficacy.xlsx"
visit_date_file = "./data/sara_testdate.xlsx"

efficacy_file_SARA = read_xlsx(efficacy_file, sheet = 1)
visit_date_SARA = read_xlsx(visit_date_file, sheet = 1)

rm(efficacy_file);rm(visit_date_file)

#### 3. Saving image file ####

save.image(file = "./data/ca_analysis.image")


#### Not used #### 
#efficacy_file_Swayneck = read_xlsx(efficacy_file, sheet = 2)
#efficacy_file_Swaywaist = read_xlsx(efficacy_file, sheet = 3)
#efficacy_file_Gaitvelocity = read_xlsx(efficacy_file, sheet = 4)
#efficacy_file_Gaitcadence = read_xlsx(efficacy_file, sheet = 5)
#efficacy_file_Gaitsupport = read_xlsx(efficacy_file, sheet = 6)
#visit_date_Wearable = read_xlsx(visit_date_file, sheet = 2)
