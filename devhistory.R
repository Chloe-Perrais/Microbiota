# Create project on Github

## Create compendium using https://github.com/benmarwick/rrtools
rrtools::use_compendium("/Users/rodenico/Documents/Pro/Articles/2022_JRLAT/Microbiota", open = FALSE)

## Add to .gitignore
usethis::use_git_ignore(".DS_Store")
usethis::use_build_ignore(".DS_Store")
## Commit and push changes
usethis::use_git(message = ":see_no_evil: Ban .DS_Store files")

## Modify DESCRIPTION file
usethis::edit_file("DESCRIPTION")
## Add information regarding the authors of the research compendium
usethis::use_git(message = ":bulb: Update documentation")

## Add license information
usethis::use_cc0_license()

## Add a read me if necesarry
## rrtools::use_readme_rmd()

dir.create("data")
dir.create("data/raw_data")
dir.create("data/derived_data")
dir.create("reports")
dir.create("plots")

## Create a R directory and a file for functions
usethis::use_r("utils-pipe")
usethis::use_r("sample_Poisson")
usethis::use_r("sample_Binomial")


## Update DESCRIPTION file
usethis::use_package("here")
usethis::use_package("ggplot2", type = "Depends")
usethis::use_package("readxl")
usethis::use_package("dplyr")
usethis::use_package("lme4")

## Update NAMESPACE file
devtools::document()

## Load all required packages
devtools::load_all()
