# 0.0 Setup ----
# 1. SnapATAC's direct dependencies, read from its GitHub DESCRIPTION file
desc <- read.dcf(url("https://raw.githubusercontent.com/r3fang/SnapATAC/master/DESCRIPTION"))
direct <- unlist(strsplit(desc[, intersect(c("Depends", "Imports", "LinkingTo"), colnames(desc))], ","))
direct <- trimws(sub("\\(.*", "", direct))
direct <- setdiff(direct[direct != ""], "R")

# 2. Add all of their dependencies, recursively
all_deps <- unique(c(direct, unlist(tools::package_dependencies(
  direct, db = installed.packages(),
  which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))

# 3. Keep only compiled packages in your personal library
mylib <- .libPaths()[1]
compiled <- all_deps[dir.exists(file.path(mylib, all_deps, "libs"))]
compiled

# 4. Reinstall them (BiocManager handles both CRAN and Bioconductor packages)
BiocManager::install(compiled, lib = mylib, force = TRUE, update = FALSE, ask = FALSE)




## 0.1 Install packages ----
.libPaths(c("~/R/R4.6.0_packages"))
options(repos = c(CRAN = "https://cloud.r-project.org")) 
install.packages("BiocManager")
BiocManager::install(c("rhdf5", "GenomicRanges", "IRanges", "S4Vectors",
                       "BiocGenerics", "edgeR", "rtracklayer"))
install.packages(c("remotes", "Matrix", "irlba", "plyr", "doSNOW", "RANN", "igraph"))
remotes::install_github("r3fang/SnapATAC")
## 0.2 Load data ----
setwd("/work/hcn4/260630_vertCons_wd/thirdTriPlacenta_snATAC/Wang_2023/paperData")
thirdTriPlacenta <- readRDS("placenta.term_atac.rds")



