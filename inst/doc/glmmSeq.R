## ----setup, include=FALSE---------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 8, fig.height = 6)
options(width=96)
library(kableExtra)

## ---- eval=FALSE------------------------------------------------------------------------------
#  install.packages("glmmSeq")

## ---- eval=FALSE------------------------------------------------------------------------------
#  devtools::install_github("KatrionaGoldmann/glmmSeq")

## ---- eval=FALSE------------------------------------------------------------------------------
#  functions = list.files("./R", full.names = TRUE)
#  invisible(lapply(functions, source))

## ---- eval=FALSE------------------------------------------------------------------------------
#  # Install CRAN packages
#  invisible(lapply(c("MASS", "car", "ggplot2", "ggpubr", "lme4", "methods",
#                     "parallel", "plotly", "stats", "gghalves"),
#                   function(p){
#                     if(! p %in% rownames(installed.packages())) {
#                       install.packages(p)
#                     }
#                     library(p, character.only=TRUE)
#                   }))
#  
#  # Install BioConductor packages
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  invisible(lapply(c("qvalue"), function(p){
#    if(! p %in% rownames(installed.packages())) BiocManager::install(p)
#    library(p, character.only=TRUE)
#  }))
#  

## ---- message=FALSE, warning=FALSE------------------------------------------------------------
library(glmmSeq)
set.seed(1234)

## ---------------------------------------------------------------------------------------------
data(PEAC_minimal_load)

## ---------------------------------------------------------------------------------------------
metadata$EULAR_binary  = NA
metadata$EULAR_binary[metadata$EULAR_6m %in%
                        c("Good responder", "Moderate responder" )] = "responder"
metadata$EULAR_binary[metadata$EULAR_6m %in% c("Non responder")] = "non_responder"
metadata = metadata[! is.na(metadata$EULAR_binary), ]

kable(head(metadata), row.names = F) %>% kable_styling()

## ---------------------------------------------------------------------------------------------
tpm = tpm[, metadata$SAMID]
kable(head(tpm)) %>% kable_styling() %>%
  scroll_box(width = "100%")

## ---------------------------------------------------------------------------------------------
disp <- apply(tpm, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
  })

head(disp)

## ---- message=FALSE---------------------------------------------------------------------------
disp  <- setNames(edgeR::estimateDisp(tpm)$tagwise.dispersion, rownames(tpm))

head(disp)

## ---- eval=FALSE------------------------------------------------------------------------------
#  dds <- DESeqDataSetFromTximport(txi = txi, colData = metadata, design = ~ 1)
#  dds <- DESeq(dds)
#  dispersions <- setNames(dispersions(dds), rownames(txi$counts))

## ---------------------------------------------------------------------------------------------
sizeFactors <- colSums(tpm)  
sizeFactors <- sizeFactors / mean(sizeFactors)  # normalise

head(sizeFactors)

## ---- eval=FALSE------------------------------------------------------------------------------
#  sizeFactors <- calcNormFactors(counts, method="TMM")

## ---- eval=FALSE------------------------------------------------------------------------------
#  sizeFactors <- estimateSizeFactorsForMatrix(counts)

## ---- warning=FALSE---------------------------------------------------------------------------
results <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
                  id = "PATID",
                  countdata = tpm,
                  metadata = metadata,
                  dispersion = disp,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  progress=TRUE,
                  cores = 1)

## ---- warning=FALSE---------------------------------------------------------------------------
results2 <- glmmSeq(~ Timepoint * EULAR_binary + (1 | PATID),
                  id = "PATID",
                  countdata = tpm,
                  metadata = metadata,
                  dispersion = disp,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  cores = 1)

## ---------------------------------------------------------------------------------------------
names(attributes(results))

## ---------------------------------------------------------------------------------------------
kable(results@modelData) %>% kable_styling()

## ---------------------------------------------------------------------------------------------
stats = data.frame(results@stats)

kable(stats[order(stats$P_Timepoint.EULAR_6m), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

## ---------------------------------------------------------------------------------------------
predict = data.frame(results@predict)
kable(predict) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

## ---------------------------------------------------------------------------------------------
results <- glmmQvals(results, pi0=1)

## ---- warning=FALSE---------------------------------------------------------------------------
MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
                     id = "PATID",
                     countdata = tpm["MS4A1", ],
                     metadata = metadata,
                     dispersion = disp,
                     verbose=FALSE)

## ---- warning=FALSE---------------------------------------------------------------------------
MS4A1fit <- glmmGene(~ Timepoint * EULAR_6m + (1 | PATID),
                     gene = "MS4A1",
                     id = "PATID",
                     countdata = tpm,
                     metadata = metadata,
                     dispersion = disp['MS4A1'])

MS4A1fit

## ---- fig.height=6, warning=FALSE-------------------------------------------------------------
plotColours <- c("skyblue", "goldenrod1", "mediumseagreen")
modColours <- c("Good responder"="dodgerblue3", 
                "Moderate responder"= "goldenrod3", 
                "Non responder"="seagreen4")
shapes <- c("Moderate responder"=19, "Good responder"= 17, "Non responder"=18)

pairedPlot(glmmResult=results,
           geneName = "IGHV3-23",
           x1Label = "Timepoint",
           x2Label="EULAR_6m",
           xTitle="Time",
           IDColumn = "PATID",
           graphics = "ggplot",
           colours = plotColours,
           shapes = shapes,
           lineColours = plotColours, 
           modelColours = modColours,
           modelLineColours = modColours,
           modelSize = 10, 
           fontSize=10,
           x2Offset = 8,
           logTransform=TRUE,
           addViolin = TRUE,
           pairedOnly = FALSE) 

## ---- fig.height=6, warning=FALSE-------------------------------------------------------------
oldpar <- par()
par(mfrow=c(1, 2))

p1 = pairedPlot(glmmResult=results2,
                geneName = "FGF14",
                x1Label = "Timepoint",
                x2Label="EULAR_binary",
                IDColumn = "PATID",
                graphics="base",
                fontSize=0.65,
                colours=c("coral", "mediumseagreen"),
                addModel=T,
                modelColours = c("coral", "mediumseagreen"),
                modelLineColours = "black",
                modelSize = 2)

p2 = pairedPlot(glmmResult=results,
                geneName = "EMILIN3",
                x1Label = "Timepoint",
                x2Label="EULAR_6m",
                IDColumn = "PATID",
                addModel=TRUE,
                graphics="base",
                fontSize=0.65,
                colours=plotColours)

par(oldpar)

## ---- message=FALSE---------------------------------------------------------------------------
library(ggpubr)

p1 <- modelPlot(results,
                "ADAM12",
                x1Label="Timepoint",
                x2Label="EULAR_6m",
                xTitle="Time",
                fontSize=8,
                x2Offset=6,
                overlap=FALSE,
                graphics="ggplot",
                colours = plotColours)

p2 <- modelPlot(results,
                "ADAM12",
                x1Label="Timepoint",
                x2Label="EULAR_6m",
                xTitle="Time",
                fontSize=8,
                x2Offset=1,
                addErrorbars = FALSE,
                overlap=TRUE,
                graphics="ggplot",
                colours = plotColours)

ggarrange(p1, p2, ncol=2, common.legend = T, legend="bottom")

## ---------------------------------------------------------------------------------------------
# Genes to label:
labels = c('MS4A1', 'FGF14', 'IL2RG', 'IGHV3-23', 'ADAM12', 'FGFRL1', 'IL36G', 
           'BLK', 'SAA1', 'CILP', 'EMILIN3', 'EMILIN2', 'IGHJ6', 
           'CXCL9', 'CXCL13')

fcPlot(glmmResult=results,
       x1Label="Timepoint",
       x2Label="EULAR_6m",
       x2Values=c("Good responder", "Non responder"),
       pCutoff=0.1,
       labels=labels,
       useAdjusted = FALSE,
       plotCutoff = 1)

## ---- fig.height=6, warning=FALSE-------------------------------------------------------------
p1<- pairedPlot(glmmResult=results,
                 geneName = "ADAM12",
                 x1Label = "Timepoint",
                 x2Label="EULAR_6m",
                 IDColumn = "PATID",
                 graphics="ggplot",
                 colours = "grey60",
                 modelColour = plotColours, 
                 modelLineColour =  plotColours, 
                 addViolins=FALSE,
                 fontSize=8,
                 logTransform=T) +
  theme(plot.subtitle=element_text(size=9))

p2 <- pairedPlot(glmmResult=results,
                 geneName = "IGHJ6",
                 x1Label = "Timepoint",
                 x2Label="EULAR_6m",
                 IDColumn = "PATID",
                 graphics="ggplot",
                 addViolins = FALSE,
                 colours = "blue",
                 fontSize=8,
                 modelSize=0.1,
                 logTransform=T) +
  theme(plot.subtitle=element_text(size=9))

ggarrange(p1, p2, ncol=2)

## ---------------------------------------------------------------------------------------------
fcPlot(glmmResult=results,
       x2Label="Timepoint",
       x1Label="EULAR_6m",
       x1Values=c("Good responder", "Non responder"),
       labels=labels,
       pCutoff=0.1,
       useAdjusted = F,
       plotCutoff = 1,
       graphics="ggplot")

## ---- fig.height=8----------------------------------------------------------------------------
maPlots <- maPlot(results,
                  x1Label="Timepoint",
                  x2Label="EULAR_6m",
                  x2Values=c("Good responder", "Non responder"),
                  colours=c('grey', 'midnightblue',
                             'mediumseagreen', 'goldenrod'),
                  labels=labels,
                  graphics="ggplot")

maPlots$combined

## ---- fig.height=8----------------------------------------------------------------------------
maPlots <- maPlot(results,
                  x2Label="Timepoint",
                  x1Label="EULAR_6m",
                  x1Values=c("Good responder", "Non responder"),
                  colours=c('grey', 'midnightblue',
                             'mediumseagreen', 'goldenrod'),
                  labels=labels,
                  graphics="ggplot")

maPlots$combined

## ---- warning=FALSE---------------------------------------------------------------------------
citation("glmmSeq")

