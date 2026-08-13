library(tidyverse)

# Load data

metadata <- readxl::read_xlsx("TCGA-TCIA/Datos/TCIA_TCGA-BRCA_09-16-2015-nbia-digest-1(1).xlsx")


# WITHOUT eliminating duplicates, obtain an object that contains only the patient ID
# the name of the protocol, the series description and the series number

metadata_protocols <- metadata[, c("Patient ID", "Protocol Name", "Series Description", "Series Number", "Image Count")] %>% 
  group_by(`Patient ID`) %>% 
  arrange(`Series Number`, .by_group = TRUE) %>% 
  as.data.frame() 
View(metadata_protocols)

# Obtain unique terms corressponding with a desired string

# Post processing terms
# Include terms that describe post processing methods 

metadata_protocols$`Series Description` %>% 
  str_subset(regex("pjn|mip|bind|cad|Processed", ignore_case = TRUE)) %>% 
  str_subset(regex("^(?!.*sub).*$", ignore_case = TRUE)) %>% 
  unique() %>% 
  cat(sep = "\n")

# Include the substraction methods: )-( to represent the observations with the format (a/b/c)-(a/d/c)

metadata_protocols$`Series Description` %>% 
  str_subset(regex("\\)-\\(", ignore_case = TRUE)) %>%
  unique() %>% 
  cat(sep = "\n")

# Other substraction methods without MIP post processing

metadata_protocols$`Series Description` %>% 
  str_subset(regex("sub", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*mip).*$", ignore_case = TRUE)) %>% 
  str_subset(regex("^(?!.*cad).*$", ignore_case = TRUE)) %>% 
  unique() %>% 
  cat(sep = "\n")

metadata_protocols$`Series Description` %>% 
  str_subset(regex("sub", ignore_case = TRUE)) %>%
  str_subset(regex("mip|cad", ignore_case = TRUE)) %>% 
  unique() %>% 
  cat(sep = "\n")




# Include de localizing images

metadata_protocols$`Series Description` %>% 
  str_subset(regex("loc", ignore_case = TRUE)) %>% 
  unique()


# T2 filtering into fat sat and not fat sat -------------------------------


# All the unique elements containing t2

all_t2 <- metadata_protocols$`Series Description` %>%
  str_subset(regex("t2|STIR|VISTA", ignore_case = TRUE)) %>%
  unique()

# T2 sat/supressed only

t2_fatsat <-  metadata_protocols$`Series Description` %>%
  str_subset(regex("t2|STIR", ignore_case = TRUE)) %>% # All T2
  str_subset(regex("sat|fs|f/s|TIRM|SPAIR|STIR|ir", ignore_case = TRUE)) %>% # Include those that are either explicit fat saturated or indirecht techiques
  str_subset(regex("^(?!.*NON FS).*$", ignore_case = TRUE)) %>% # Ignore those explicitly stated as Non fat saturated
  str_subset(regex("^(?!.*T SAG FSE).*$", ignore_case = FALSE)) %>% # Ignore those that oinclude FSE as a whole since that is a technique that is not fat saturated
  str_subset(regex("^(?!.*T sag fse).*$", ignore_case = FALSE)) %>% # Same as before
  str_subset(regex("^(?!.*T2 fse ax bilat).*$", ignore_case = FALSE)) %>% # Same as before
  unique() # The length of this object is 12

t2_nonfs <- metadata_protocols$`Series Description` %>%
  str_subset(regex("t2|NON FS", ignore_case = TRUE)) %>%  # Obtain t2
  str_subset(regex("^(?!.*sat).*$", ignore_case = TRUE)) %>%  # This and the next 4 lines ignore either explicit fat saturated descriptions or a technique of fat saturation
  str_subset(regex("^(?!.*f/s).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*TIRM).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*SPAIR).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*\\bir\\b).*$", ignore_case = TRUE)) %>% # This one is to ignore ir method of fat saturation
  str_subset(regex("^(?!.*(?<!NON )\\bfs\\b)(?!.*f/s)(?!.*sat)(?!.*TIRM)(?!.*SPAIR).*$", ignore_case = TRUE)) %>% # This is so taht it does not eliminate the NON FS object
  unique() # The length of this object is 14 which summed with previous comletes the 26 uniques of t2

# Undetermined T Weight 3D

und_t2 <- metadata_protocols$`Series Description` %>%
  str_subset(regex("vista|space|cube", ignore_case = TRUE)) %>% # T2 VISTA
  str_subset(regex("^(?!.*CAD).*$")) %>%
  unique()



# T1 filtering into pre and post contrast -------------------------------


# All the unique elements containing t1

metadata_protocols$`Series Description` %>% 
  str_subset(regex("t1", ignore_case = TRUE)) %>% 
  unique() %>%
length() # Length is of 28




# T1 pre contrast fat saturated

t1_fs <- metadata_protocols$`Series Description` %>% 
  str_subset(regex("t1", ignore_case = TRUE)) %>% # All T1
  str_subset(regex("fat sat|fs |spair", ignore_case = TRUE)) %>% # Include those registered as fat saturated
  str_subset(regex("^(?!.*post).*$", ignore_case = TRUE)) %>% # Eliminate those registered as post contrast
  str_subset(regex("^(?!.*non fat).*$", ignore_case = TRUE)) %>% # eliminate those that explicitly state non fat
  str_subset(regex("^(?!.*not fat).*$", ignore_case = TRUE)) %>% # eliminate those that explicitly state non fat
  unique() # Length of 3
t1_fs %>% 
  cat(sep = "\n")

# T1 pre contrast not fat saturated

metadata_protocols$`Series Description` %>% 
  str_subset(regex("t1", ignore_case = TRUE)) %>% # All T1
  str_subset(regex("^(?!.*post).*$", ignore_case = TRUE)) %>% # Include those that are not post contrast
  str_subset(regex(paste0("^(?!.*", t1_fs[1], ").*$"), ignore_case = TRUE)) %>% # This and the next line are to eliminate the terms already defiened as fat sat
  str_subset(regex(paste0("^(?!.*", t1_fs[2], ").*$"), ignore_case = TRUE)) %>% # done this way to not filter out accidentaly the "non fat" objects
  str_subset(regex(paste0("^(?!.*", t1_fs[3], ").*$"), ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*ffe).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*Flash).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*SPAIR).*$", ignore_case = TRUE)) %>%
  unique() %>% # Length of 16
  cat(sep = "\n")




# DCE (Dynamic Contrast-Enhanced) filtering -------------------------------

# Clean DCE sequences only (filtering out the post-processing and technical artifacts)

dce_expl_term <- metadata_protocols$`Series Description` %>%
  str_subset(regex("vibrant|thrive|dynamic|flash|Views|post|ffe|\\d+\\s*min", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*loc).*$", ignore_case = TRUE)) %>%        # Eliminate localizers
  str_subset(regex("^(?!.*mip).*$", ignore_case = TRUE)) %>%        # Eliminate Maximum Intensity Projections
  str_subset(regex("^(?!.*cad).*$", ignore_case = TRUE)) %>%        # Eliminate Computer Aided Detection
  str_subset(regex("^(?!.*collection).*$", ignore_case = TRUE)) %>% # Eliminate collections
  str_subset(regex("^(?!.*bind).*$", ignore_case = TRUE)) %>%       # Eliminate bind data
  str_subset(regex("^(?!.*scout).*$", ignore_case = TRUE)) %>%      # Eliminate scout views
  str_subset(regex("^(?!.*pjn).*$", ignore_case = TRUE)) %>%        # Eliminate projection lines
  str_subset(regex("^(?!.*\\)-\\().*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*sub).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*US).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*INPHASE).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*FAT:).*$")) %>%
  str_subset(regex("^(?!.*MOCO).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*interviews).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*pilot).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*hires).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*radial).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*post-contrast).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*water).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*sag).*$", ignore_case = TRUE)) %>%

  unique()

dce_expl_term2 <- metadata_protocols$`Series Description` %>%
  str_subset(regex("\\)\\s*-\\s*\\(", ignore_case = TRUE))

dce_expl_term %>%
  cat(sep = "\n")

dce_expl_pt <- metadata_protocols %>%
  filter((`Series Description` %in% dce_expl_term))


dce_expl_pt_compmetadata <- metadata_protocols %>%
  filter((`Patient ID` %in% dce_expl_pt$`Patient ID`))

length(unique(dce_expl_pt$`Patient ID`))

dce_prob <- metadata_protocols %>%
  filter(!(`Patient ID` %in% dce_expl_pt$`Patient ID`),
         grepl("\\)-\\(", `Series Description`))

unique(dce_prob$`Patient ID`)

dce_prob_pt <- metadata_protocols %>%
  filter((`Patient ID` %in% dce_prob$`Patient ID`))


dce_t2 <- metadata_protocols %>%
  filter((`Patient ID` %in% dce_expl_pt$`Patient ID`) &
           (`Series Description` %in% c(t2_fatsat, t2_nonfs, und_t2)))

z <- unique(dce_t2$`Patient ID`)


for (i in unique(dce_prob$`Patient ID`)) {
  cat("\n")
  
  print(i)
  
  dce_prob_pt_loop <- dce_prob_pt %>%
    filter(`Patient ID` == as.character(i))
  
  
  
  cat(dce_prob_pt_loop$`Series Number`[ grepl("Bind\\(", dce_prob_pt_loop$`Series Description`) == TRUE], sep = "\n")
  
  
  
  cat(dce_prob_pt_loop$`Series Description`[ grepl("Bind\\(", dce_prob_pt_loop$`Series Description`) == TRUE], sep = "\n")
  
  cat(dce_prob_pt_loop$`Series Description`[ grepl("\\)-\\(", dce_prob_pt_loop$`Series Description`) == TRUE], sep = "\n")
  cat("\n")
  
}




# DWI ---------------------------------------------------------------------

dwi <- metadata_protocols$`Series Description` %>% 
  str_subset(regex("DWI|diff", ignore_case = TRUE)) %>% # All T1
  str_subset(regex("^(?!.*cad).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*Apparent).*$", ignore_case = TRUE)) %>%
  str_subset(regex("^(?!.*adc).*$", ignore_case = TRUE)) %>%
  unique() 

metadata_protocols$`Series Description` %>% 
  str_subset(regex("adc|diff_ORTHOGONALXYZ_ADC|Diffusion", ignore_case = TRUE)) %>% # All T1
  str_subset(regex("^(?!.*cad).*$", ignore_case = TRUE)) %>%
  unique() %>%  # The length of this object is 7
  cat(sep = "\n")
