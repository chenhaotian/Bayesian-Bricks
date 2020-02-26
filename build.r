library(roxygen2)
library(devtools)
library(usethis)

## remove previous build
system("rm -r man")
system("rm NAMESPACE DESCRIPTION")

## create skeleton files
create_package("../bbricks",rstudio = FALSE)

## generate man/*.Rd files and NAMESPACE
devtools::document()

## add to .Rbuildignore to ignore non-package files
usethis::use_build_ignore(c("LICENSE.md",".travis.yml","build.r",".#build.r","README.raw.md","notes_pictures"))

## Edit DESCRIPTIONN
## Version Reference:
##    http://r-pkgs.had.co.nz/release.html
##    Version format:  major.minor.patch.dev
## add to Collate: when there's new files added
Description <- paste0("Package: bbricks
Title: bbricks
Version: 0.1.0.9000
Authors@R: 
    person(given = \"Haotian\",
           family = \"Chen\",
           role = c(\"aut\", \"cre\"),
           email = \"chenhaotian.jtt@gmail.com\",
           comment = structure(\"https://orcid.org/0000-0001-9751-2093\", .Names = \"ORCID\"))
Description: Basic building blocks in Bayesian network modeling.
License: What license it uses
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

install.packages("../bbricks_0.1.0.9000.tar.gz") #different from install()
## install()                                        #different from install.packages()

library(bbricks)

remove.packages("bbricks")
## uninstall()
## remove.packages("bbricks")
