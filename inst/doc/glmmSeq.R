## ----setup, include=FALSE-------------------------------------------------------------
options(width=88)
library(kableExtra)

## ---- eval=FALSE----------------------------------------------------------------------
#  install.packages("glmmSeq")

## ---- eval=FALSE----------------------------------------------------------------------
#  devtools::install_github("KatrionaGoldmann/glmmSeq")

## ---- eval=FALSE----------------------------------------------------------------------
#  functions = list.files("./R", full.names = TRUE)
#  invisible(lapply(functions, source))

## ---- eval=FALSE----------------------------------------------------------------------
#  # Install CRAN packages
#  invisible(lapply(c("MASS", "car", "ggplot2", "ggpubr", "lme4","lmerTest",
#                     "methods", "parallel", "plotly", "pbapply", "pbmcapply"),
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

## ---- message=FALSE, warning=FALSE----------------------------------------------------
library(glmmSeq)
set.seed(1234)

## -------------------------------------------------------------------------------------
data(PEAC_minimal_load)

## -------------------------------------------------------------------------------------
metadata$EULAR_binary  = NA
metadata$EULAR_binary[metadata$EULAR_6m %in%
                        c("Good", "Moderate" )] = "responder"
metadata$EULAR_binary[metadata$EULAR_6m %in% c("Non-response")] = "non_responder"

kable(head(metadata), row.names = F) %>% kable_styling()

## -------------------------------------------------------------------------------------

kable(head(tpm)) %>% kable_styling() %>%
  scroll_box(width = "100%")

## -------------------------------------------------------------------------------------
disp <- apply(tpm, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
  })

head(disp)

## ---- eval=FALSE----------------------------------------------------------------------
#  disp  <- setNames(edgeR::estimateDisp(tpm)$tagwise.dispersion, rownames(tpm))
#  
#  head(disp)

## ---- eval=FALSE----------------------------------------------------------------------
#  dds <- DESeqDataSetFromTximport(txi = txi, colData = metadata, design = ~ 1)
#  dds <- DESeq(dds)
#  dispersions <- setNames(dispersions(dds), rownames(txi$counts))

## -------------------------------------------------------------------------------------
sizeFactors <- colSums(tpm)  
sizeFactors <- sizeFactors / mean(sizeFactors)  # normalise to mean = 1

head(sizeFactors)

## ---- eval=FALSE----------------------------------------------------------------------
#  sizeFactors <- calcNormFactors(counts, method="TMM")

## ---- eval=FALSE----------------------------------------------------------------------
#  sizeFactors <- estimateSizeFactorsForMatrix(counts)

## ---- warning=FALSE-------------------------------------------------------------------
results <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
                   countdata = tpm,
                   metadata = metadata,
                   dispersion = disp,
                   progress = TRUE)

## ---- warning=FALSE-------------------------------------------------------------------
results2 <- glmmSeq(~ Timepoint * EULAR_binary + (1 | PATID),
                    countdata = tpm,
                    metadata = metadata,
                    dispersion = disp)

## -------------------------------------------------------------------------------------
names(attributes(results))

## -------------------------------------------------------------------------------------
kable(results@modelData) %>% kable_styling()

## -------------------------------------------------------------------------------------
stats <- summary(results)

kable(stats[order(stats[, 'P_Timepoint:EULAR_6m']), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

## -------------------------------------------------------------------------------------
summary(results, gene = "MS4A1")

## -------------------------------------------------------------------------------------
predict = data.frame(results@predict)
kable(predict) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

## -------------------------------------------------------------------------------------
results <- glmmQvals(results)

## ---- warning=FALSE-------------------------------------------------------------------
logtpm <- log2(tpm + 1)
lmmres <- lmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
                   maindata = logtpm,
                   metadata = metadata,
                   progress = TRUE)
summary(lmmres, "MS4A1")

## ---- warning=FALSE-------------------------------------------------------------------
MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
                     countdata = tpm["MS4A1", , drop = FALSE],
                     metadata = metadata,
                     dispersion = disp,
                     verbose = FALSE)

## ---- warning=FALSE-------------------------------------------------------------------
fit <- glmmRefit(results, gene = "MS4A1")
fit

## ---- warning=FALSE-------------------------------------------------------------------
library(emmeans)

emmeans(fit, ~ Timepoint | EULAR_6m)
emmip(fit, ~ Timepoint | EULAR_6m)

## ---- fig.height=6, warning=FALSE-----------------------------------------------------
plotColours <- c("skyblue", "goldenrod1", "mediumseagreen")
modColours <- c("dodgerblue3", "goldenrod3", "seagreen4")
shapes <- c(17, 19, 18)

ggmodelPlot(results,
            geneName = "IGHV3-23",
            x1var = "Timepoint",
            x2var="EULAR_6m",
            xlab="Time",
            colours = plotColours,
            shapes = shapes,
            lineColours = plotColours, 
            modelColours = modColours,
            modelSize = 10)

## ---- fig.height=6, warning=FALSE-----------------------------------------------------
oldpar <- par(mfrow=c(1, 2))

modelPlot(results2,
          geneName = "FGF14",
          x1var = "Timepoint",
          x2var="EULAR_binary",
          fontSize=0.65,
          colours=c("coral", "mediumseagreen"),
          modelColours = c("coral", "mediumseagreen"),
          modelLineColours = "black",
          modelSize = 2)

modelPlot(results,
          geneName = "EMILIN3",
          x1var = "Timepoint",
          x2var = "EULAR_6m",
          colours = plotColours,
          plab = c("time", "response", "time:response"),
          addModel = FALSE)

par(oldpar)

## ---- message=FALSE-------------------------------------------------------------------
library(ggpubr)

p1 <- ggmodelPlot(results,
                  "ADAM12",
                  x1var="Timepoint",
                  x2var="EULAR_6m",
                  xlab="Time",
                  addPoints = FALSE,
                  colours = plotColours)

p2 <- ggmodelPlot(results,
                  "EMILIN3",
                  x1var="Timepoint",
                  x2var="EULAR_6m",
                  xlab="Time",
                  fontSize=8,
                  x2Offset=1,
                  addPoints = FALSE,
                  colours = plotColours)

ggarrange(p1, p2, ncol=2, common.legend = T, legend="bottom")

## ---- eval = FALSE--------------------------------------------------------------------
#  r4ra_glmm <- glmmSeq(~ time * drug + (1 | Patient_ID),
#                         countdata = tpmdata, metadata,
#                         dispersion = dispersions, cores = 8, removeSingles = T)
#  r4ra_glmm <- glmmQvals(r4ra_glmm)
#  labels = c(..)  # Genes to label
#  fcPlot(r4ra_glmm, x1var = "time", x2var = "drug", graphics = "plotly",
#         pCutoff = 0.05, useAdjusted = TRUE,
#         labels = labels,
#         colours = c('grey', 'green3', 'gold3', 'blue'))

## ----fcplot, echo = FALSE, message=FALSE, fig.align='center', out.width='80%', out.extra='style="border: 0;"'----
knitr::include_graphics("r4ra_glmm_fcplot.png")

## ---- fig.height=8--------------------------------------------------------------------
labels = c('MS4A1', 'FGF14', 'IL2RG', 'IGHV3-23', 'ADAM12', 'IL36G', 
           'BLK', 'SAA1', 'CILP', 'EMILIN3', 'EMILIN2', 'IGHJ6', 
           'CXCL9', 'CXCL13')
maPlots <- maPlot(results,
                  x1var="Timepoint",
                  x2var="EULAR_6m",
                  x2Values=c("Good", "Non-response"),
                  colours=c('grey', 'midnightblue',
                             'mediumseagreen', 'goldenrod'),
                  labels=labels,
                  graphics="ggplot")

maPlots$combined

## ---- warning=FALSE-------------------------------------------------------------------
citation("glmmSeq")

