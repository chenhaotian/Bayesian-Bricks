library(roxygen2)
library(devtools)
library(usethis)

VERSION <- "0.1.0.9000"
VERSION <- "0.1.0"

## remove previous build
system("rm -r man")
system("rm -r vignettes Meta doc")
system("rm NAMESPACE DESCRIPTION .Rbuildignore")

## create skeleton files
create_package("../bbricks",rstudio = FALSE)

## generate man/*.Rd files and NAMESPACE
devtools::document()

## add to .Rbuildignore to ignore non-package files
## "NEWS.md" is supported by CRAN, no need to exclude
usethis::use_build_ignore(c("LICENSE.md","README.md",".travis.yml","build.r",".#build.r","README.raw.md","notes_pictures","cran-comments.md"))
## usethis::use_build_ignore(c("README.md",".travis.yml","build.r",".#build.r","README.raw.md","cran-comments.md"))

## Edit DESCRIPTIONN
## Version Reference:
##    http://r-pkgs.had.co.nz/release.html
##    Version format:  major.minor.patch.dev
## add to Collate: when there's new files added
Description <- paste0("Type: Package
Package: bbricks
Title: Bayesian methods and graphical model structures for statistical modeling
Version: ",VERSION,"
Authors@R: 
    person(given = \"Haotian\",
           family = \"Chen\",
           role = c(\"aut\", \"cre\"),
           email = \"chenhaotian.jtt@gmail.com\",
           comment = structure(\"0000-0001-9751-2093\", .Names = \"ORCID\"))
Description: A class of frequently used Bayesian parametric and nonparametric model structures and the corresponding inference tasks. Structures include Gaussian and Normal-Inverse-Wishart conjugate structure, Gaussian and Normal-Inverse-Gamma conjugate structure, Categorical and Dirichlet conjugate structure, Dirichlet Process on positive integers, Dirichlet Process in general, Hierarchical Dirichlet Process and so on; Tasks include updating posteriors, calculating marginal likelihood, calculating posterior predictives, calculating MAP estimates and so on.
License: What license it uses
URL: https://github.com/chenhaotian/Bayesian-Bricks
BugReports: https://github.com/chenhaotian/Bayesian-Bricks/issues
Encoding: UTF-8
LazyData: true
RoxygenNote: ",packageVersion("roxygen2"),"
Collate: 
    'Bayesian_Bricks.r'
    'Categorical_Inference.r'
    'Gaussian_Inference.r'
    'Dirichlet_Process.r'
    'MCMC.r'
    'bbricks-package.R'
    'testData.r'
")
write(Description,file = "DESCRIPTION")

## Add license, run this will update the corresponding line in DESCRIPTION automatically
use_mit_license(name="Haotian Chen")

## Config travis CI
travis <- "language: r
r:
  - oldrel
  - release
  - devel
  - bioc-devel
  - bioc-release
r_packages:
  - stats
"
write(travis,file = ".travis.yml")


## add vignette
## reference: https://cran.r-project.org/web/packages/R.rsp/vignettes/R_packages-RSP_vignettes.html
use_vignette("bbricks-getting-started")

readme <- readLines("README.raw.md")
pics <- gsub("\\).*$","",gsub("^.*\\(","",grep("\\./notes_pictures/",readme,value = TRUE)))
if(length(pics)>0){
    if(length(dir("vignettes/notes_pictures/"))==0)
        system("mkdir vignettes/notes_pictures")
    system(paste("cp",paste(pics,collapse = " "),"vignettes/notes_pictures/"))
}
picsNotInline <- grep("^\\!.*\\./notes_pictures/.*\\)$",readme)
readme[picsNotInline] <- paste0(readme[picsNotInline],"{width=100%}")

readme <- gsub("```R","```{r,eval=FALSE}",readme)

writeLines(c("---
title: \"bbricks-getting-started\"
output: rmarkdown::html_vignette
vignette: >
  %\\VignetteIndexEntry{bbricks: Getting Started}
  %\\VignetteEngine{knitr::rmarkdown}
  %\\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = \"#>\"
)
```
",readme),"vignettes/bbricks-getting-started.Rmd")

system("rm .#* *~")
build_vignettes()

## very important step, must make sure there's no  warnings and errors
## remove stupic emacs tmp file before check()
system("rm .#* *~")
devtools::check()



## somehow the S3 methods are not exported... export them mannually
## have to be after check(), otherwise it will be overwritten
namespace <- readLines("./NAMESPACE")
S3classes <- grep("^S3method\\(",namespace,value = TRUE)
S3classes <- unique(gsub(",.*$","",gsub("S3method\\(","",S3classes)))
funs <- gsub(".rd$","",dir("./man"),ignore.case = TRUE)
funs <- na.omit(sapply(funs,function(f){
    if(any(sapply(S3classes,function(s3){grepl(paste0("^",s3,"\\."),f)}))) f
    else as.character(NA)
},USE.NAMES = FALSE))
funs <- paste0("export(",funs,")")
namespace <- c(namespace,setdiff(funs,namespace)) #don't use union()
writeLines(namespace,"./NAMESPACE")


## build the package
## remove stupid emacs tmp file before build()
system("rm .#* *~")
devtools::build()

## devtools::release()

install.packages(paste0("../bbricks_",VERSION,".tar.gz")) #different from install()
## install()                                        #different from install.packages()

library(bbricks)

remove.packages("bbricks")
## uninstall()
## remove.packages("bbricks")


## check from rhub
check_rhub()

## check from win builder
check_win_devel()
check_win_release()
check_win_oldrelease()

## If you are ready to sumit it to CRAN, run relsease, release will do following things:
## • Confirm that the package passes ‘R CMD check’ on relevant
##    platforms
## • Confirm that important files are up-to-date
## • Build the package
## • Submit the package to CRAN, using comments in
##    "cran-comments.md"
release_questions <- function() {
  c("Have you set the correct version number?",
    "Have you removed the irrelevant code blocks?",
    "Have you add {width=100%} to each inluded image?",
    "Have you add all R files to DESCRIPTION?",
    "Have you removed the unfinished lines from vignette?")
}

devtools::release()

