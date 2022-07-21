library(rbioclim)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rgbif)
library(dplyr)
library(moiR)
################################################################################
# Set map limits
xlim = c(-200, +200)
ylim = c(0, +90)
# Download present data w2
now<-rbioclim::getData(name="worldclim2",var='bio', res=2.5, path = "~/")
now<- now %>% cropenvironment(.,xlim,ylim)
# Download future data Max Planck model
fut<-rbioclim::getData(name="CMIP5",var='bio', res=2.5,model='MP', year=50, rcp=85,path = "~/")
fut<- fut %>% cropenvironment(.,xlim,ylim)
# Fix names in future dataset
names(fut) <- names(now)
# Fix units of temp in future dataset (wordlclim2 temp is not Cx10)
for(i in 1:11) fut[[i]]<-fut[[i]]/10
################################################################################
# Arabidopsis data
# GBIF.org (16 June 2022) GBIF Occurrence Download  https://doi.org/10.15468/dl.pda86b
# https://api.gbif.org/v1/occurrence/download/request/0357669-210914110416597.zip
# Download the data manually, too large for GBIF
gbif_<-read.delim("0357669-210914110416597.csv",header = T, sep='\t', fill = T)
colnames(gbif_)
gbif<-gbif_ %>%
  dplyr::filter(occurrenceStatus=='PRESENT',
                decimalLongitude > -15, 
                decimalLatitude > 15)
# Subsample 1 sample per degree square
gbifsub<-gbif %>%
  dplyr::mutate(latround=round(decimalLatitude),
                lonround=round(decimalLongitude)
  ) %>%
  dplyr::group_by(latround,lonround) %>%
  sample_n(1)
coords<- gbif[,c("decimalLongitude","decimalLatitude")]
################################################################################
# Quick plot
# ggplot(gbif)+
ggplot(gbifsub)+
  geom_point(aes(y=decimalLatitude, x=decimalLongitude))
################################################################################
# Extract variables
bioclim<-raster::extract(x = now ,coords)
bioclimfut<-raster::extract(x = fut ,coords)
hist(bioclimfut[,1] - bioclim[,1])
hist(bioclimfut[,5] - bioclim[,5])
################################################################################
# Niche model
library(dismo)
library(rJava)
mod<-maxent(x = now, p=coords)
sdm<-predict(mod, now) 
plot(sdm)
sdmfut<-predict(mod, fut) 
# visualize
plot(sdmfut)
plot((sdmfut>0.5)-(sdm>0.5))
plot(sdmfut-sdm)
change=crop(sdmfut-sdm,extent(c(-15, 53),c(30,65)))
pdf(file = "arabidopsis-change.pdf")
plot(change,col=brewer.pal(11,"RdBu"))
dev.off()
hist(sdm-sdmfut)
# Save info
save(file = "mod.rda",mod)
save(file = "sdm.rda",sdm)
save(file = "sdmfut.rda",sdmfut)
