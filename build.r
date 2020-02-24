library(roxygen2)
library(devtools)
library(usethis)

system("rm -r man")
system("rm NAMESPACE DESCRIPTION")


create_package("../bbricks",rstudio = FALSE)

document()

Rbuildignore <- "^LICENSE\\.md$
^\\.travis\\.yml$
^build\\.r$
^\\.#build\\.r$
"
write(Rbuildignore,file = ".Rbuildignore")

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
    'bbricks-package.R'
    'testData.r'
")

write(Description,file = "DESCRIPTION")

use_mit_license(name="Haotian Chen")

## devtools::use_travis()                  #add to .Rbuildignore

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

system("rm .#* *~")
devtools::build()

system("rm .#* *~")

## very important step
## must make sure there's no  warnings and errors
devtools::check()

## devtools::release()

install()

library(bbricks)

uninstall()
## remove.packages("bbricks")
