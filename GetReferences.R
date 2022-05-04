# Citations for Maize Genetic Background Microbiome Paper
# Corey Schultz - 2022
library(tidyr)
library(purrr)
library(stringr)

#### get all the package names:

# Get all files
getwd()
WD = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts"
setwd(WD)
files <- list.files(path=WD, pattern=".R", full.names=T, recursive = TRUE)
print(files)

# make new folder 
new.folder <- "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/DuplicateText"

dir.create(new.folder)

setwd(new.folder)

file.copy(files,new.folder)

dup_files <- list.files(path = new.folder, pattern = ".R", all.files = T)

new_r_scripts <- sub(pattern="\\.R$", replacement=".txt", x=dup_files)
file.rename(from = dup_files, to = new_r_scripts)

# Read in Scripts
scripts <- lapply(new_r_scripts, function(x)readChar(x, file.info(x)$size))

#turn the list given from the above function into a vector
script_vector <- unlist(scripts)

# Regex
packages <- sort(unique(unlist(noquote(stringr::str_extract_all(script_vector, "library\\(.*")))))
packages <- stringr::str_extract_all(packages, "(?<=\\().+?(?=\\))")
packages <- sort(str_replace_all(packages, "[[:punct:]]", ""))

# sub any packages with interior punctiuation here. 
packages <- gsub("datatable","data.table",packages)

print(noquote(packages))

# Mass citation
capture.output(noquote(packages) %>%
  map(citation) %>%
  print(style = "text"), file = "package_references.txt")

# append R info!
write(capture.output(citation()), file = "package_references.txt", append = T)
write(version$version.string, file = "package_references.txt", append = T)

