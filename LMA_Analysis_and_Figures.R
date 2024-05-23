########Streamlined tests of LMA stuff

#Load in Packages
library(tidyverse)
library(plotly)
library(philentropy)
library(graphics)
library(ggraph)
library(dendextend)
library(viridis)
library(vegan)
library(cowplot)
#library(ggbiome)
library(whitIt)
library(plotbiomes)
library(ggrepel)
library(GGally)
library(RColorBrewer)
library(patchwork)
library(multcomp)
library(FD)
library(gridGraphics)
library(ggpp)
library(sp)

#Add functions
colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,hdContext$...1)
    treatment=hdContext[ligne,2];
    if(treatment=="Modern"){col_treatment="blue"};if(treatment=="Fossil"){col_treatment="red"}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,col=col_treatment,lab.font=1,lab.cex=.5))
  }
  return(n)
}

regressionp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an lm")
  f <- summary(modelobject$fstatistic)
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  return(p)
}

#Read in datasheets
fossilSpecimens <- read_csv("fossilSpecimens.csv") #Individual fossil leaf data, most has LMA info and depositional environment
fossilHistograms <- read.csv("fossilHistograms.csv", header = FALSE, row.names = 1)
fossilSpeciesOnly <- read_csv("fossilSpeciesOnly.csv") %>% dplyr::select(c(1:4)) %>%
  mutate(logLMA = log10(MeanLMA))
modernSpecimens <- readxl::read_xlsx("modernSpecimens.xlsx", sheet = "Data Compilation") %>% slice(-1)
siteCoords <- read.csv("siteCoords.csv")
modernSiteClimate <- read_csv("modernSiteClimate.csv")
modernSpecimenTaxonomy <- read_csv("modernSpecimenTaxonomy.csv")
fossilSiteClimate <- read_csv("fossilSiteClimate.csv")
fossilMorphology <- read_csv("fossilMorphology.csv")
fossilMorphology2 <- read_csv("fossilMorphology2.csv", skip = 1)

#Clean and filter modern specimen data and calculate site means
filteredModernLMA <- modernSpecimens %>%
  dplyr::select(1:9, 14, 32:40, 116) %>%
  filter(!is.na(LMA)) %>%
  filter(`Angio/Gymno` == "dicot") %>%
  filter(`Plant life form` %in% c("Tree", "Erect dwarf shrub", "Prostrate dwarf shrub", "Liana", "Small Tree", "Shrub", "Low to high shrub")) %>%
  filter(`Leaf type` %in% c(NA, "Broad")) %>%
  group_by(Database, Site, Lat, Lon, Alt, Species, ApprovedSpecies, SummaryEvDec) %>%
  summarise(MeanLMA = mean(LMA), logLMA = log10(MeanLMA), n = n()) %>%
  group_by(Database, Site, Lat, Lon) %>%
  mutate("numbSpecies" = n()) %>%
  ungroup() %>%
  filter(numbSpecies>=10) %>%
  dplyr::select(!numbSpecies)
filteredModernSites <- filteredModernLMA %>%
  group_by(Database, Site, Lat, Lon, Alt) %>%
  summarise(MeanLMA = mean(MeanLMA), logLMA = mean(logLMA), "numSpecies" = n(), decNum = replace_na(length(which(SummaryEvDec=="Deciduous")),0), evNum = replace_na(length(which(SummaryEvDec=="Evergreen")),0), naNum = replace_na(length(which(SummaryEvDec=="NA")),0)) %>%
  filter(numSpecies >= 10) %>%
  ungroup()

######Part 1: Maps
###########1A: Map of Fossil localities
fossilMap <- ggplot() + geom_map(data = map_data("world"), map = map_data("world"), aes(long, lat, map_id = region), 
                                 fill = "antiquewhite", color = "black") + 
  geom_point(data = filter(siteCoords, Set == "Fossil"), 
             aes(as.numeric(Lon), as.numeric(Lat)), size = 2, fill = "red", shape = 21) +
  labs(title = "Fossil Localities") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))

###########1B: Map of Modern day Localities
modernMap <- ggplot() + geom_map(data = map_data("world"), map = map_data("world"), aes(long, lat, map_id = region), 
                                 fill = "antiquewhite", color = "black") + 
  geom_point(data = filter(siteCoords, Set == "Modern"), 
             aes(as.numeric(Lon), as.numeric(Lat)), size = 2, fill = "blue", shape = 21) +
  labs(title = "Modern Localities") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))

#Figure 1
figure1 <- plot_grid(fossilMap, modernMap, ncol = 1, labels = c("A", "B"))
rm(fossilMap, modernMap)

######Part 2: Whittaker Diagrams
###########2A: Whittaker Diagram of Fossil localities
fossilWhit <- ggplot() +  geom_polygon(data = Whittaker_biomes, aes(x = temp_c, y = precp_cm*10, fill = biome),
                                       color = "gray98", size = 1) + theme_bw() +
  scale_fill_manual(name   = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = c(Ricklefs_colors)) +
  labs(x = expression("Mean Annual Temperature ("*degree*C*")"), y = "Mean Annual Precipitation (mm)", title = "Fossil Localities") +
  geom_point(data = fossilSiteClimate, aes(x = MAT, y = MAP * 10), size = 2, stroke = .3, alpha = .75, fill = "red", shape = 21) + 
  theme(legend.position = "none")
legend <- get_legend(fossilWhit + theme(legend.position = "right"))

###########2B: Whittaker Diagram of Modern localities
whitClim <- filteredModernSites[1:5] %>%
  left_join(modernSiteClimate, by = c("Database", "Site", "Lat", "Lon", "Alt"))

modernWhit <- ggplot() +  geom_polygon(data = Whittaker_biomes, aes(x = temp_c, y = precp_cm * 10, fill = biome),
                                       color = "gray98", size = 1) + theme_bw() +
  scale_fill_manual(name   = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = c(Ricklefs_colors)) +
  labs(x = expression("Mean Annual Temperature ("*degree*C*")"), y = "Mean Annual Precipitation (mm)", title = "Modern Localities") +
  geom_point(data = whitClim, aes(x = newMAT, y = newMAP), size = 2, stroke = .3, alpha = .75, fill = "blue", shape = 21) + 
  theme(legend.position = "none")

#Figure 2
figure2 <- plot_grid(fossilWhit, legend, modernWhit, ncol = 2, rel_widths = c(2,1), labels = c("A", "", "B"))
rm(fossilWhite, legend, modernWhit)


######Part 3: LMA Regressions
###########3A: First build a cleaned fossil dataset
filteredModernLMA <- filteredModernLMA %>%
  left_join(whitClim, by = c("Database", "Site", "Lat", "Lon", "Alt"))
filteredModernSites <- filteredModernSites %>%
  left_join(whitClim, by = c("Database", "Site", "Lat", "Lon", "Alt"))

temp <- fossilSpecimens %>%
  group_by(Site, Family, Genus, Species, ID) %>%
  summarize("MeanLMA" = 10^(log10(mean(PW)^2/mean(LA)) * 0.382 + 3.07), "samples" = n()) %>%
  mutate("logLMA" = log10(MeanLMA))
tempSites <- temp %>%
  group_by(Site) %>%
  summarise("n" = n())
depEnv <- fossilSpecimens %>%
  group_by(Site, DepEnv) %>%
  summarise("k" = n())
filteredFossilLMA <- left_join(temp, tempSites, by = "Site") %>%
  bind_rows(fossilSpeciesOnly) %>%
  left_join(fossilSiteClimate, by = c("Site", "n")) %>%
  left_join(whitClim, by = "Site") %>%
  group_by(Site) %>%
  mutate("numbSpecies" = n()) %>%
  ungroup() %>%
  filter(numbSpecies >= 10) %>%
  dplyr::select(!numbSpecies)
filteredFossilSites <- filteredFossilLMA %>%
  group_by(Site) %>%
  summarize(MeanLMA = mean(MeanLMA), logLMA = mean(logLMA),"numSpecies" = n()) %>%
  filter(numSpecies >=10) %>%
  ungroup() %>%
  left_join(fossilSiteClimate, by = "Site") %>%
  mutate(MAP = MAP * 10)
rm(temp, depEnv, tempSites)

##Constrain modern dataset to same MAT and MAP range as fossil dataset
lowRangeModern <- filteredModernSites %>%
  filter(newMAT > 7.29) %>%
  filter(newMAT < 28.01) %>%
  filter(newMAP > 589) %>%
  filter(newMAP < 3241)

#Test quadratic correlations between site-mean LMA and key climate variables
filteredModernLMAquad <- filteredModernSites %>%
  dplyr::select(c(Site, MeanLMA, logLMA, TempSeasonality, gridVPD, gridRSDS, newMAT, newMAP)) %>%
  mutate(temp2 = TempSeasonality^2, vpd2 = gridVPD^2, rsds2 = gridRSDS^2, mat2 = newMAT ^2, map2 = log10(newMAP)^2, log = logLMA^2)

quadSeason <- ggplot(data = filteredModernSites, aes(x = logLMA, y = TempSeasonality)) + geom_point(color = "blue") + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) + geom_smooth(method = "lm", formula = y ~ x, color = "red") +  ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "red"),                                                                                                                                            data = data.frame(x = 0.05, y = 0.1, label = "R^2 = 0.32")) +
  ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "blue"), 
                      data = data.frame(x = 0.05, y = 0.05, label = "R^2 = 0.32")) + theme_bw() +
  labs(y = expression("Temperature Seasonality (SD * 100)"), x = expression("log10(Leaf Mass per Area (g/m"^2*"))"))
quadVPD <- ggplot(data = filteredModernSites, aes(x = logLMA, y = gridVPD)) + geom_point(color = "blue") + geom_smooth(method = "lm", formula = y ~ x + I(x^2))+ geom_smooth(method = "lm", formula = y ~ x, color = "red") +  ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "red"), data = data.frame(x = 0.05, y = 0.1, label = "R^2 = 0.13")) + ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "blue"), data = data.frame(x = 0.05, y = 0.05, label = "R^2 = 0.12")) + theme_bw() +
  labs(y = "Vapor Pressure Deficit (Pa)", x = expression("log10(Leaf Mass per Area (g/m"^2*"))"))
quadRSDS <- ggplot(data = filteredModernSites, aes(x = logLMA, y = gridRSDS)) + geom_point(color = "blue") + geom_smooth(method = "lm", formula = y ~ x + I(x^2))+ geom_smooth(method = "lm", formula = y ~ x, color = "red") +  ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "red"), data = data.frame(x = 0.05, y = 0.1, label = "R^2 = 0.18")) + ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "blue"), data = data.frame(x = 0.05, y = 0.05, label = "R^2 = 0.19")) + theme_bw() +
  labs(y = expression("RSDS (MJ/m"^2*"d)"), y = expression("log10(Leaf Mass per Area (g/m"^2*"))"))
quadNAM <- ggplot(data = filteredModernSites, aes(x = logLMA, y = newMAT)) + geom_point(color = "blue") + geom_smooth(method = "lm", formula = y ~ x + I(x^2))+ geom_smooth(method = "lm", formula = y ~ x, color = "red") +  ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "red"), data = data.frame(x = 0.05, y = 0.1, label = "R^2 = 0.18")) + ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "blue"), data = data.frame(x = 0.05, y = 0.05, label = "R^2 = 0.18")) + theme_bw() +
  labs(y = expression("Mean Annual Temperature ("*degree*C*")"), x = expression("log10(Leaf Mass per Area (g/m"^2*"))"))
quadMAP <- ggplot(data = filteredModernSites, aes(x = logLMA, y = log10(newMAP))) + geom_point(color = "blue") + geom_smooth(method = "lm", formula = y ~ x + I(x^2))+ geom_smooth(method = "lm", formula = y ~ x, color = "red") +  ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "red"), data = data.frame(x = 0.05, y = 0.1, label = "R^2 = 0.04")) + ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label, color = "blue"), data = data.frame(x = 0.05, y = 0.05, label = "R^2 = 0.06")) + theme_bw() +
  labs(y = "Mean Annual Precipitation (mm)", x = expression("log10(Leaf Mass per Area (g/m"^2*"))"))
#Appendix S9
appendixS9 <- quadNAM + quadMAP + plot_spacer() + quadRSDS + quadSeason + quadVPD + 
  plot_annotation(tag_levels = 'A')
rm(quadNAM, quadMAP, quadRSDS, quadVPD, quadSeason)

#Regression models for the linear and quadratic relationships
siteMeanRegressions <- list(
  "MAT Linear" = summary(lm(logLMA ~newMAT, data = filteredModernSites)),
  "MAT Quadratic" = summary(lm(logLMA ~newMAT + mat2, data = filteredModernLMAquad)),
  "MAP Linear" = summary(lm(logLMA ~log10(newMAP), data = filteredModernSites)),
  "MAP Quadratic" = summary(lm(logLMA ~log10(newMAP) + map2, data = filteredModernLMAquad)),
  "RSDS Linear" = summary(lm(logLMA ~gridRSDS, data = filteredModernSites)),
  "RSDS Quadratic" = summary(lm(logLMA ~gridRSDS + rsds2, data = filteredModernLMAquad)),
  "Temp Seasonality Linear" = summary(lm(logLMA ~TempSeasonality, data = filteredModernSites)),
  "Temp Seasonality Quadratic" = summary(lm(logLMA ~TempSeasonality + temp2, data = filteredModernLMAquad)),
  "VPD Linear" = summary(lm(logLMA ~gridVPD, data = filteredModernSites)),
  "VPD Quadratic" = summary(lm(logLMA ~gridVPD + vpd2, data = filteredModernLMAquad))
)

###########3B: LMA vs MAT
matPlot <- ggplot() + geom_point(data = filteredModernLMA, aes(x = newMAT, y = MeanLMA), alpha = 0.1, size = .9) + 
  geom_point(data = filteredModernSites, aes(x = newMAT, y = MeanLMA), color = "blue") +
  geom_smooth(data = filteredModernSites, aes(x = newMAT, y = MeanLMA), method = "lm", color = "blue", se = FALSE) +
  geom_point(data = filteredFossilLMA, aes(x = MAT, y = MeanLMA), alpha = .1, size = .9) +
  geom_point(data = filteredFossilSites, aes(x = MAT, y = MeanLMA), color = "red") +
  geom_smooth(data = filteredFossilSites, aes(x = MAT, y = MeanLMA), method = "lm", color = "red", se = FALSE) +
  theme_bw() + scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  labs(x = expression("Mean Annual Temperature ("*degree*C*")"), y = expression("Leaf Mass per Area (g/m"^2*")"))
###########3C: LMA vs MAP
mapPlot <- ggplot() + geom_point(data = filteredModernLMA, aes(x = newMAP, y = MeanLMA), alpha = 0.1, size = .9) + 
  geom_point(data = filteredFossilLMA, aes(x = MAP * 10, y = MeanLMA), alpha = .1, size = .9) +
  geom_point(data = filteredModernSites, aes(x = newMAP, y = MeanLMA), color = "blue") +
  geom_smooth(data = filteredModernSites, aes(x = newMAP, y = MeanLMA), method = "lm", color = "blue", se = FALSE) +
  geom_point(data = filteredFossilSites, aes(x = MAP, y = MeanLMA), color = "red") +
  geom_smooth(data = filteredFossilSites, aes(x = MAP, y = MeanLMA), method = "lm", color = "red", se = FALSE) +
  theme_bw() + scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  #scale_x_continuous(trans = 'log10') +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()) +
  labs(x = "Mean Annual Precipitation (mm)", y = expression("Leaf Mass per Area (g/m"^2*")"))
###########3D: LMA vs gridVPD
vpdPlot <- ggplot() + geom_point(data = filteredModernLMA, aes(x = gridVPD, y = MeanLMA), alpha = 0.1, size = .9) + 
  geom_point(data = filteredModernSites, aes(x = gridVPD, y = MeanLMA), color = "blue") +
  geom_smooth(data = filteredModernSites, aes(x = gridVPD, y = MeanLMA), method = "lm", color = "blue", se = FALSE) +
  theme_bw() + scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()) +
  labs(x = "Vapor Pressure Deficit (Pa)", y = expression("Leaf Mass per Area (g/m"^2*")"))
###########3E: LMA vs gridPET
petPlot <- ggplot() + geom_point(data = filteredModernLMA, aes(x = gridPET, y = MeanLMA), alpha = 0.1, size = .9) + 
  geom_point(data = filteredModernSites, aes(x = gridPET, y = MeanLMA), color = "blue") +
  geom_smooth(data = filteredModernSites, aes(x = gridPET, y = MeanLMA), method = "lm", color = "blue", se = FALSE) +
  theme_bw() + scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()) +
  labs(x = expression("Potential Evapotranpsiration (kg/m"^2*")"), y = expression("Leaf Mass per Area (g/m"^2*")"))
###########3F: LMA vs surface downwelling shortwave radiation
rsdsPlot <- ggplot() + geom_point(data = filteredModernLMA, aes(x = gridRSDS, y = MeanLMA), alpha = 0.1, size = .9) + 
  geom_point(data = filteredModernSites, aes(x = gridRSDS, y = MeanLMA), color = "blue") +
  geom_smooth(data = filteredModernSites, aes(x = gridRSDS, y = MeanLMA), method = "lm", color = "blue", se = FALSE) +
  theme_bw() + scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  labs(x = expression("Surface Downwelling Shortwave Radiation (MJ/m"^2*"d)"), y = expression("Leaf Mass per Area (g/m"^2*")"))
############3G: LMA vs CO2 (Fossil)
co2Plot <- ggplot() + geom_point(data = filteredFossilLMA, aes(x = co2, y = MeanLMA), alpha = 0.1, size = .9) +
  geom_point(data = filteredFossilSites, aes(x = co2, y = MeanLMA), color = "red",) +
  geom_point(data = filteredModernLMA, aes(x = 1000, y = MeanLMA), alpha = 0) +
  geom_smooth(data = filteredFossilSites, aes(x = co2, y = MeanLMA), color = "red", method = "lm", se = FALSE) +
  theme_bw() + scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()) +
  labs(x = expression("Atmospheric CO"[2]*" (ppm)"), y = expression("Leaf Mass per Area (g/m"^2*")"))
###########3H: LMA vs TempSeasonality
seasonPlot <- ggplot() + geom_point(data = filteredModernLMA, aes(x = TempSeasonality, y = MeanLMA), alpha = 0.1, size = .9) + 
  geom_point(data = filteredModernSites, aes(x = TempSeasonality, y = MeanLMA), color = "blue") +
  geom_smooth(data = filteredModernSites, aes(x = TempSeasonality, y = MeanLMA), method = "lm", color = "blue", se = FALSE) +
  theme_bw() + scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()) +
  labs(x = expression("Temperature Seasonality (SD * 100)"), y = expression("Leaf Mass per Area (g/m"^2*")"))
climateVars <- filteredModernSites %>%
  dplyr::select(c(7,12:21, 23:29, 32:43, 46:47)) %>%
  mutate(newMAP = log10(newMAP), PrecipWettestMonth = log10(PrecipWettestMonth), PrecipDriestMonth = log10(PrecipDriestMonth + 1), PrecipWettestQuarter = log10(PrecipWettestQuarter), PrecipDriestQuarter = log10(PrecipDriestQuarter + 1), PrecipWarmestQuarter = log10(PrecipWarmestQuarter), PrecipColdestQuarter = log10(PrecipColdestQuarter))
lowRangeClimateVars <- lowRangeModern %>%
  dplyr::select(c(7,12:21, 23:29, 32:43, 46:47)) %>%
  mutate(newMAP = log10(newMAP), PrecipWettestMonth = log10(PrecipWettestMonth), PrecipDriestMonth = log10(PrecipDriestMonth + 1), PrecipWettestQuarter = log10(PrecipWettestQuarter), PrecipDriestQuarter = log10(PrecipDriestQuarter + 1), PrecipWarmestQuarter = log10(PrecipWarmestQuarter), PrecipColdestQuarter = log10(PrecipColdestQuarter))
#Figure3
figure3 <- matPlot + mapPlot + co2Plot + rsdsPlot + seasonPlot + vpdPlot + plot_layout(nrow = 2) + 
  plot_annotation(tag_levels = 'A')
rm(vpdPlot, seasonPlot, matPlot, mapPlot, co2Plot, rsdsPlot, petPlot)

#Appendix S8 - Covariance between climate variables
appendixS8 <- cor(climateVars, use = "pairwise.complete.obs")

#######Get coefficients for site linear regressions
tf <- data.frame("R" = 0, "p" = 0, "Slope" = 0, "Variable" = "", "t" = 0)
for (i in 1) {
  tt <- summary(lm(logLMA ~gridCMI, data = filteredModernSites))
  tf[1,1] <- sqrt(tt$r.squared)
  tf[1,2] <- tt$coefficients[2,4]
  tf[1,3] <- tt$coefficients[2,1]
  tf[1,4] <- rownames(tt$coefficients)[2]
  tf[1,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridNPP, data = filteredModernSites))
  tf[2,1] <- sqrt(tt$r.squared)
  tf[2,2] <- tt$coefficients[2,4]
  tf[2,3] <- tt$coefficients[2,1]
  tf[2,4] <- rownames(tt$coefficients)[2]
  tf[2,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridPET, data = filteredModernSites))
  tf[3,1] <- sqrt(tt$r.squared)
  tf[3,2] <- tt$coefficients[2,4]
  tf[3,3] <- tt$coefficients[2,1]
  tf[3,4] <- rownames(tt$coefficients)[2]
  tf[3,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridVPD, data = filteredModernSites))
  tf[4,1] <- sqrt(tt$r.squared)
  tf[4,2] <- tt$coefficients[2,4]
  tf[4,3] <- tt$coefficients[2,1]
  tf[4,4] <- rownames(tt$coefficients)[2]
  tf[4,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridFCF, data = filteredModernSites))
  tf[5,1] <- sqrt(tt$r.squared)
  tf[5,2] <- tt$coefficients[2,4]
  tf[5,3] <- tt$coefficients[2,1]
  tf[5,4] <- rownames(tt$coefficients)[2]
  tf[5,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridRSDS, data = filteredModernSites))
  tf[6,1] <- sqrt(tt$r.squared)
  tf[6,2] <- tt$coefficients[2,4]
  tf[6,3] <- tt$coefficients[2,1]
  tf[6,4] <- rownames(tt$coefficients)[2]
  tf[6,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridHURS, data = filteredModernSites))
  tf[7,1] <- sqrt(tt$r.squared)
  tf[7,2] <- tt$coefficients[2,4]
  tf[7,3] <- tt$coefficients[2,1]
  tf[7,4] <- rownames(tt$coefficients)[2]
  tf[7,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridCLT, data = filteredModernSites))
  tf[8,1] <- sqrt(tt$r.squared)
  tf[8,2] <- tt$coefficients[2,4]
  tf[8,3] <- tt$coefficients[2,1]
  tf[8,4] <- rownames(tt$coefficients)[2]
  tf[8,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridSCD, data = filteredModernSites))
  tf[9,1] <- sqrt(tt$r.squared)
  tf[9,2] <- tt$coefficients[2,4]
  tf[9,3] <- tt$coefficients[2,1]
  tf[9,4] <- rownames(tt$coefficients)[2]
  tf[9,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridGSL, data = filteredModernSites))
  tf[10,1] <- sqrt(tt$r.squared)
  tf[10,2] <- tt$coefficients[2,4]
  tf[10,3] <- tt$coefficients[2,1]
  tf[10,4] <- rownames(tt$coefficients)[2]
  tf[10,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridSFCWind, data = filteredModernSites))
  tf[11,1] <- sqrt(tt$r.squared)
  tf[11,2] <- tt$coefficients[2,4]
  tf[11,3] <- tt$coefficients[2,1]
  tf[11,4] <- rownames(tt$coefficients)[2]
  tf[11,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridGST, data = filteredModernSites))
  tf[12,1] <- sqrt(tt$r.squared)
  tf[12,2] <- tt$coefficients[2,4]
  tf[12,3] <- tt$coefficients[2,1]
  tf[12,4] <- rownames(tt$coefficients)[2]
  tf[12,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~DiurnalRange, data = filteredModernSites))
  tf[13,1] <- sqrt(tt$r.squared)
  tf[13,2] <- tt$coefficients[2,4]
  tf[13,3] <- tt$coefficients[2,1]
  tf[13,4] <- rownames(tt$coefficients)[2]
  tf[13,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~Isothermality, data = filteredModernSites))
  tf[14,1] <- sqrt(tt$r.squared)
  tf[14,2] <- tt$coefficients[2,4]
  tf[14,3] <- tt$coefficients[2,1]
  tf[14,4] <- rownames(tt$coefficients)[2]
  tf[14,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~TempSeasonality, data = filteredModernSites))
  tf[15,1] <- sqrt(tt$r.squared)
  tf[15,2] <- tt$coefficients[2,4]
  tf[15,3] <- tt$coefficients[2,1]
  tf[15,4] <- rownames(tt$coefficients)[2]
  tf[15,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MaxTempWarmestMonth, data = filteredModernSites))
  tf[16,1] <- sqrt(tt$r.squared)
  tf[16,2] <- tt$coefficients[2,4]
  tf[16,3] <- tt$coefficients[2,1]
  tf[16,4] <- rownames(tt$coefficients)[2]
  tf[16,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MinTempColdestMonth, data = filteredModernSites))
  tf[17,1] <- sqrt(tt$r.squared)
  tf[17,2] <- tt$coefficients[2,4]
  tf[17,3] <- tt$coefficients[2,1]
  tf[17,4] <- rownames(tt$coefficients)[2]
  tf[17,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~TempAnnualRange, data = filteredModernSites))
  tf[18,1] <- sqrt(tt$r.squared)
  tf[18,2] <- tt$coefficients[2,4]
  tf[18,3] <- tt$coefficients[2,1]
  tf[18,4] <- rownames(tt$coefficients)[2]
  tf[18,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MeanTempWettestQuarter, data = filteredModernSites))
  tf[19,1] <- sqrt(tt$r.squared)
  tf[19,2] <- tt$coefficients[2,4]
  tf[19,3] <- tt$coefficients[2,1]
  tf[19,4] <- rownames(tt$coefficients)[2]
  tf[19,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MeanTempDriestQuarter, data = filteredModernSites))
  tf[20,1] <- sqrt(tt$r.squared)
  tf[20,2] <- tt$coefficients[2,4]
  tf[20,3] <- tt$coefficients[2,1]
  tf[20,4] <- rownames(tt$coefficients)[2]
  tf[20,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MeanTempColdestQuarter, data = filteredModernSites))
  tf[21,1] <- sqrt(tt$r.squared)
  tf[21,2] <- tt$coefficients[2,4]
  tf[21,3] <- tt$coefficients[2,1]
  tf[21,4] <- rownames(tt$coefficients)[2]
  tf[21,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipWettestMonth), data = filteredModernSites))
  tf[22,1] <- sqrt(tt$r.squared)
  tf[22,2] <- tt$coefficients[2,4]
  tf[22,3] <- tt$coefficients[2,1]
  tf[22,4] <- rownames(tt$coefficients)[2]
  tf[22,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipDriestMonth +1), data = filteredModernSites))
  tf[23,1] <- sqrt(tt$r.squared)
  tf[23,2] <- tt$coefficients[2,4]
  tf[23,3] <- tt$coefficients[2,1]
  tf[23,4] <- rownames(tt$coefficients)[2]
  tf[23,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~PrecipSeasonality, data = filteredModernSites))
  tf[24,1] <- sqrt(tt$r.squared)
  tf[24,2] <- tt$coefficients[2,4]
  tf[24,3] <- tt$coefficients[2,1]
  tf[24,4] <- rownames(tt$coefficients)[2]
  tf[24,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipWettestQuarter), data = filteredModernSites))
  tf[25,1] <- sqrt(tt$r.squared)
  tf[25,2] <- tt$coefficients[2,4]
  tf[25,3] <- tt$coefficients[2,1]
  tf[25,4] <- rownames(tt$coefficients)[2]
  tf[25,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipDriestQuarter+1), data = filteredModernSites))
  tf[26,1] <- sqrt(tt$r.squared)
  tf[26,2] <- tt$coefficients[2,4]
  tf[26,3] <- tt$coefficients[2,1]
  tf[26,4] <- rownames(tt$coefficients)[2]
  tf[26,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipWarmestQuarter), data = filteredModernSites))
  tf[27,1] <- sqrt(tt$r.squared)
  tf[27,2] <- tt$coefficients[2,4]
  tf[27,3] <- tt$coefficients[2,1]
  tf[27,4] <- rownames(tt$coefficients)[2]
  tf[27,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipColdestQuarter), data = filteredModernSites))
  tf[28,1] <- sqrt(tt$r.squared)
  tf[28,2] <- tt$coefficients[2,4]
  tf[28,3] <- tt$coefficients[2,1]
  tf[28,4] <- rownames(tt$coefficients)[2]
  tf[28,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~newMAT, data = filteredModernSites))
  tf[29,1] <- sqrt(tt$r.squared)
  tf[29,2] <- tt$coefficients[2,4]
  tf[29,3] <- tt$coefficients[2,1]
  tf[29,4] <- rownames(tt$coefficients)[2]
  tf[29,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(newMAP), data = filteredModernSites))
  tf[30,1] <- sqrt(tt$r.squared)
  tf[30,2] <- tt$coefficients[2,4]
  tf[30,3] <- tt$coefficients[2,1]
  tf[30,4] <- rownames(tt$coefficients)[2]
  tf[30,5] <- tt$coefficients[2,3]
}
#Appendix S7
siteRegressions <- tf %>% dplyr::select(Variable, Slope, t, R, p) %>%
  mutate(sig = (p<=0.05))

#######Get coefficients for species linear regressions
tf <- data.frame("R" = 0, "p" = 0, "Slope" = 0, "Variable" = "", "t" = 0)
for (i in 1) {
  tt <- summary(lm(logLMA ~gridCMI, data = filteredModernLMA))
  tf[1,1] <- sqrt(tt$r.squared)
  tf[1,2] <- tt$coefficients[2,4]
  tf[1,3] <- tt$coefficients[2,1]
  tf[1,4] <- rownames(tt$coefficients)[2]
  tf[1,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridNPP, data = filteredModernLMA))
  tf[2,1] <- sqrt(tt$r.squared)
  tf[2,2] <- tt$coefficients[2,4]
  tf[2,3] <- tt$coefficients[2,1]
  tf[2,4] <- rownames(tt$coefficients)[2]
  tf[2,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridPET, data = filteredModernLMA))
  tf[3,1] <- sqrt(tt$r.squared)
  tf[3,2] <- tt$coefficients[2,4]
  tf[3,3] <- tt$coefficients[2,1]
  tf[3,4] <- rownames(tt$coefficients)[2]
  tf[3,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridVPD, data = filteredModernLMA))
  tf[4,1] <- sqrt(tt$r.squared)
  tf[4,2] <- tt$coefficients[2,4]
  tf[4,3] <- tt$coefficients[2,1]
  tf[4,4] <- rownames(tt$coefficients)[2]
  tf[4,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridFCF, data = filteredModernLMA))
  tf[5,1] <- sqrt(tt$r.squared)
  tf[5,2] <- tt$coefficients[2,4]
  tf[5,3] <- tt$coefficients[2,1]
  tf[5,4] <- rownames(tt$coefficients)[2]
  tf[5,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridRSDS, data = filteredModernLMA))
  tf[6,1] <- sqrt(tt$r.squared)
  tf[6,2] <- tt$coefficients[2,4]
  tf[6,3] <- tt$coefficients[2,1]
  tf[6,4] <- rownames(tt$coefficients)[2]
  tf[6,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridHURS, data = filteredModernLMA))
  tf[7,1] <- sqrt(tt$r.squared)
  tf[7,2] <- tt$coefficients[2,4]
  tf[7,3] <- tt$coefficients[2,1]
  tf[7,4] <- rownames(tt$coefficients)[2]
  tf[7,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridCLT, data = filteredModernLMA))
  tf[8,1] <- sqrt(tt$r.squared)
  tf[8,2] <- tt$coefficients[2,4]
  tf[8,3] <- tt$coefficients[2,1]
  tf[8,4] <- rownames(tt$coefficients)[2]
  tf[8,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridSCD, data = filteredModernLMA))
  tf[9,1] <- sqrt(tt$r.squared)
  tf[9,2] <- tt$coefficients[2,4]
  tf[9,3] <- tt$coefficients[2,1]
  tf[9,4] <- rownames(tt$coefficients)[2]
  tf[9,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridGSL, data = filteredModernLMA))
  tf[10,1] <- sqrt(tt$r.squared)
  tf[10,2] <- tt$coefficients[2,4]
  tf[10,3] <- tt$coefficients[2,1]
  tf[10,4] <- rownames(tt$coefficients)[2]
  tf[10,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridSFCWind, data = filteredModernLMA))
  tf[11,1] <- sqrt(tt$r.squared)
  tf[11,2] <- tt$coefficients[2,4]
  tf[11,3] <- tt$coefficients[2,1]
  tf[11,4] <- rownames(tt$coefficients)[2]
  tf[11,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~gridGST, data = filteredModernLMA))
  tf[12,1] <- sqrt(tt$r.squared)
  tf[12,2] <- tt$coefficients[2,4]
  tf[12,3] <- tt$coefficients[2,1]
  tf[12,4] <- rownames(tt$coefficients)[2]
  tf[12,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~DiurnalRange, data = filteredModernLMA))
  tf[13,1] <- sqrt(tt$r.squared)
  tf[13,2] <- tt$coefficients[2,4]
  tf[13,3] <- tt$coefficients[2,1]
  tf[13,4] <- rownames(tt$coefficients)[2]
  tf[13,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~Isothermality, data = filteredModernLMA))
  tf[14,1] <- sqrt(tt$r.squared)
  tf[14,2] <- tt$coefficients[2,4]
  tf[14,3] <- tt$coefficients[2,1]
  tf[14,4] <- rownames(tt$coefficients)[2]
  tf[14,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~TempSeasonality, data = filteredModernLMA))
  tf[15,1] <- sqrt(tt$r.squared)
  tf[15,2] <- tt$coefficients[2,4]
  tf[15,3] <- tt$coefficients[2,1]
  tf[15,4] <- rownames(tt$coefficients)[2]
  tf[15,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MaxTempWarmestMonth, data = filteredModernLMA))
  tf[16,1] <- sqrt(tt$r.squared)
  tf[16,2] <- tt$coefficients[2,4]
  tf[16,3] <- tt$coefficients[2,1]
  tf[16,4] <- rownames(tt$coefficients)[2]
  tf[16,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MinTempColdestMonth, data = filteredModernLMA))
  tf[17,1] <- sqrt(tt$r.squared)
  tf[17,2] <- tt$coefficients[2,4]
  tf[17,3] <- tt$coefficients[2,1]
  tf[17,4] <- rownames(tt$coefficients)[2]
  tf[17,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~TempAnnualRange, data = filteredModernLMA))
  tf[18,1] <- sqrt(tt$r.squared)
  tf[18,2] <- tt$coefficients[2,4]
  tf[18,3] <- tt$coefficients[2,1]
  tf[18,4] <- rownames(tt$coefficients)[2]
  tf[18,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MeanTempWettestQuarter, data = filteredModernLMA))
  tf[19,1] <- sqrt(tt$r.squared)
  tf[19,2] <- tt$coefficients[2,4]
  tf[19,3] <- tt$coefficients[2,1]
  tf[19,4] <- rownames(tt$coefficients)[2]
  tf[19,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MeanTempDriestQuarter, data = filteredModernLMA))
  tf[20,1] <- sqrt(tt$r.squared)
  tf[20,2] <- tt$coefficients[2,4]
  tf[20,3] <- tt$coefficients[2,1]
  tf[20,4] <- rownames(tt$coefficients)[2]
  tf[20,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~MeanTempColdestQuarter, data = filteredModernLMA))
  tf[21,1] <- sqrt(tt$r.squared)
  tf[21,2] <- tt$coefficients[2,4]
  tf[21,3] <- tt$coefficients[2,1]
  tf[21,4] <- rownames(tt$coefficients)[2]
  tf[21,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipWettestMonth), data = filteredModernLMA))
  tf[22,1] <- sqrt(tt$r.squared)
  tf[22,2] <- tt$coefficients[2,4]
  tf[22,3] <- tt$coefficients[2,1]
  tf[22,4] <- rownames(tt$coefficients)[2]
  tf[22,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipDriestMonth + 1), data = filteredModernLMA))
  tf[23,1] <- sqrt(tt$r.squared)
  tf[23,2] <- tt$coefficients[2,4]
  tf[23,3] <- tt$coefficients[2,1]
  tf[23,4] <- rownames(tt$coefficients)[2]
  tf[23,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~PrecipSeasonality, data = filteredModernLMA))
  tf[24,1] <- sqrt(tt$r.squared)
  tf[24,2] <- tt$coefficients[2,4]
  tf[24,3] <- tt$coefficients[2,1]
  tf[24,4] <- rownames(tt$coefficients)[2]
  tf[24,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipWettestQuarter), data = filteredModernLMA))
  tf[25,1] <- sqrt(tt$r.squared)
  tf[25,2] <- tt$coefficients[2,4]
  tf[25,3] <- tt$coefficients[2,1]
  tf[25,4] <- rownames(tt$coefficients)[2]
  tf[25,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipDriestQuarter+1), data = filteredModernLMA))
  tf[26,1] <- sqrt(tt$r.squared)
  tf[26,2] <- tt$coefficients[2,4]
  tf[26,3] <- tt$coefficients[2,1]
  tf[26,4] <- rownames(tt$coefficients)[2]
  tf[26,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipWarmestQuarter), data = filteredModernLMA))
  tf[27,1] <- sqrt(tt$r.squared)
  tf[27,2] <- tt$coefficients[2,4]
  tf[27,3] <- tt$coefficients[2,1]
  tf[27,4] <- rownames(tt$coefficients)[2]
  tf[27,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(PrecipColdestQuarter), data = filteredModernLMA))
  tf[28,1] <- sqrt(tt$r.squared)
  tf[28,2] <- tt$coefficients[2,4]
  tf[28,3] <- tt$coefficients[2,1]
  tf[28,4] <- rownames(tt$coefficients)[2]
  tf[28,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~newMAT, data = filteredModernLMA))
  tf[29,1] <- sqrt(tt$r.squared)
  tf[29,2] <- tt$coefficients[2,4]
  tf[29,3] <- tt$coefficients[2,1]
  tf[29,4] <- rownames(tt$coefficients)[2]
  tf[29,5] <- tt$coefficients[2,3]
  tt <- summary(lm(logLMA ~log10(newMAP), data = filteredModernLMA))
  tf[30,1] <- sqrt(tt$r.squared)
  tf[30,2] <- tt$coefficients[2,4]
  tf[30,3] <- tt$coefficients[2,1]
  tf[30,4] <- rownames(tt$coefficients)[2]
  tf[30,5] <- tt$coefficients[2,3]
}
speciesRegressions <- tf %>% dplyr::select(Variable, Slope, t, R, p) %>%
  mutate(sig = (p<=0.05)) 

appendixS7 <- list("Species LMA Regressions" = speciesRegressions, "Site LMA Regressions" = siteRegressions)
rm(tf, tt, speciesRegressions, siteRegressions)

######Part 4: Introduce Distributions, clustering methods, show clusters, and average distribution curves for each cluster.
#Make probability distributions for modern sites
histStarter <- filteredModernLMA %>%
  ungroup() %>%
  left_join(dplyr::select(filteredModernSites, c("Database", "Site", "Lat", "Lon", "Alt", "numSpecies")), by = c("Database", "Site", "Lat", "Lon", "Alt")) %>%
  filter(!is.na(numSpecies)) %>%
  dplyr::select(Site, logLMA) %>%
  pivot_wider(values_from = logLMA, names_from = Site) %>%
  t()
#PDFs with gaussian kernal between min and max LMA value of .65 and 2.9
histStarter2 <- unlist(density(as.numeric(unlist(histStarter[1,])), kernel = "gaussian", from = .65, to = 2.9, na.rm = TRUE)[2])
for (i in 1:length(histStarter)) {
  histStarter2 <- bind_rows(histStarter2, unlist(density(as.numeric(unlist(histStarter[i,])), kernel = "gaussian", from = .65, to = 2.9, na.rm = TRUE)[2]))
}
histStarter2 <- histStarter2 %>%
  slice(-1)
#Bind fossil site distributions to table
for (i in 1:nrow(fossilHistograms)) {
  histStarter2 <- bind_rows(histStarter2, unlist(density(as.numeric(log10(fossilHistograms[i,])), kernel = "gaussian", from = .65, to = 2.9, na.rm = TRUE)[2]))
}
rows <- c(rownames(histStarter), rownames(fossilHistograms))
#Form distance matrix
dm <- JSD(as.matrix(histStarter2))
rownames(dm) <- rows
hd <- hclust(as.dist(dm), method = "complete")
hdContext <- bind_cols(hd$labels, NULL)
hdContext[1:length(histStarter),2] <- "Modern"
hdContext[(length(histStarter) + 1):311 ,2] <- "Fossil"

dend <- as.dendrogram(hd)
#Make dendrogram based on LMA distributions
dL <- dendrapply(dend, colLab)
#Divide dendrogram into 7 clusters
logClusters <- cutree(dL,7)
tempFrame <- data.frame(logClusters = c(1,2,3,4,5,6,7), clusterName = c("4. Middle Right", "5. Low Right", "2. Low Left", "7. Middle Rightmost", "3. Middle Left", "6. Upper Right", "1. Middle Leftmost"))
logClusters <- rownames_to_column(as.data.frame(logClusters)) %>%
  left_join(tempFrame, by = "logClusters")
colnames(logClusters)[1] <- "Site"
filteredFossilLMA <- filteredFossilLMA %>%
  left_join(logClusters, by = "Site")
filteredModernLMA <- filteredModernLMA %>%
  left_join(logClusters, by = "Site")
filteredFossilSites <- filteredFossilSites %>%
  left_join(logClusters, by = "Site")
filteredModernSites <- filteredModernSites %>%
  left_join(logClusters, by = "Site")

###########4A: Cluster cladogram

palette <-  c('#e377c2',  # raspberry yogurt pink/ middle leftmost
              '#1f77b4',  # muted blue/low left
              '#2ca02c',  # cooked asparagus green/ middle left
              '#9467bd',  # muted purple/ middle right
              '#d62728',  # brick red/ low right
              '#8c564b',  # chestnut brown/ upper right
              '#ff7f0e'   # safety orange/ middle rightmost
              
)
par(mar = c(1,4,4,2))
dL %>%
  set("labels_cex", 0.3) %>%
  set("branches_k_color", k = 7, value = c(    '#e377c2',  # raspberry yogurt pink/middleleftmost
                                               '#8c564b',  # chestnut brown/upper right
                                               '#ff7f0e',  # safety orange/middle rightmost
                                               '#9467bd',  # muted purple/middle right
                                               '#d62728',   # 7brick red/ low right
                                               '#1f77b4',  # muted blue/low left
                                               '#2ca02c'  # 6cooked asparagus green/middle left
  )) %>% 
  plot(main = "Fossil and Modern Distributions", leaflab = "none")
legend("topright", 
       legend = c("Fossil Site" , "Modern Site"), 
       col = c("red", "blue"), 
       pch = c(20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0.05, 0.05))
cladogram <- recordPlot()

rm(dend, hd, hdContext, histStarter, dm, rows, i, tempFrame)

###########4B: Average curves for each cluster
tempFossils <- filteredFossilLMA %>% ungroup() %>%
  dplyr::select(Site, Species = ID, logLMA, MeanLMA, logClusters, clusterName) %>%
  mutate(Set = "Fossil")
allSpecies <- filteredModernLMA %>% ungroup() %>%
  dplyr::select(Site, Species = ApprovedSpecies, logLMA, MeanLMA, logClusters, clusterName) %>%
  mutate(Set = "Modern") %>%
  bind_rows(tempFossils) %>%
  filter(logClusters != 0)
rm(tempFossils)
distPlot <- ggplot() + 
  geom_line(data = allSpecies, aes(MeanLMA, group = Site, color = clusterName), alpha = 0.3, stat = "density") +
  geom_line(data = allSpecies, aes(MeanLMA, color = clusterName), size = 2, stat = "density") +
  scale_color_manual(values = palette) + scale_x_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  labs(x = expression("Leaf Mass per Area (g/m"^2*")"), y = "Density") +
  theme_bw() + theme(legend.position = "none")
#Figure 4
figure4 <- plot_grid(cladogram, distPlot, nrow = 2, ncol = 1, labels = c("A", "B"))
rm(cladogram, distPlot, dL)

######Part 5: Clusters with key types of data.
###########5A: Clusters with the family NMDS - not particularly affected by the taxonomic composition of the flora!
familyCircle <- modernSpecimenTaxonomy %>% 
  filter(!is.na(family)) %>%
  dplyr::select(Database, Site, family, n) %>%
  pivot_wider(names_from = family, values_from = n, names_sort = TRUE, values_fn = mean) %>%
  mutate_at(c(3:length(.)), ~replace_na(.,0)) %>%
  dplyr::select(!Database) %>%
  group_by(Site) %>%
  summarize(across(.cols = everything(), .fns = sum)) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  column_to_rownames(var = "Site")


fc2 <- familyCircle %>%
  mutate(sum = rowSums(across())) %>%
  filter(sum>14) %>%
  dplyr::select(!sum)

temp <- summarize(fc2, across(everything(), ~sum(., is.na(.), 0)))

fc2 <- fc2 %>%
  bind_rows(temp) %>%
  select_if(function(.) last(.) >1) %>%
  head(-1)

familyMatch <- fc2 %>%
  rownames_to_column(var = "Site") %>%
  left_join(whitClim, by = "Site") %>%
  left_join(logClusters, by = "Site") %>%
  dplyr::select(Site, logClusters, clusterName)

taxOrdPlot <- metaMDS(fc2, k =2, distance = "raup", trymax =  150)
taxonomyOrdinationScores <- as.data.frame(scores(taxOrdPlot)$sites) %>%
  rownames_to_column(var = "Site") %>%
  left_join(familyMatch, by = "Site") %>%
  dplyr::select(Site, NMDS1, NMDS2, clusterName)
#Figure 7 - Taxonomy NMDS
par(mar = c(5,4,4,2))
ordiplot(taxOrdPlot, choices = c(1,2), display = "sites")
ordispider(taxOrdPlot, choices =c(1,2), groups = familyMatch$clusterName, label = FALSE, col = c('#1f77b4',  # muted blue/low left
                                                                                                '#2ca02c',  # cooked asp/middle left
                                                                                                '#9467bd',  #purple/middle right
                                                                                                '#d62728'  # brick red/ low right
), show.groups = c("2. Low Left", "3. Middle Left", "5. Low Right", "4. Middle Right"))
legend("bottomleft", col = c('#1f77b4',  # muted blue/low left
                             '#d62728',  # brick red/ low right
                             '#2ca02c',  # cooked asp/middle left
                             '#9467bd'), lty =1, bty = "n", legend = c("Low Left", "Low Right", "Middle Left", "Middle Right"))

figure7 <- recordPlot()

#Calculate distance between clusters using ANOSIM
attach(familyMatch)
taxOrdPlot.ano <- anosim(fc2, grouping = clusterName, distance = "raup")
#plot(taxOrdPlot.ano)
#summary(taxOrdPlot.ano)
detach(familyMatch)
rm(temp, fc2, familyMatch)

###########5B: Clusters with Whittaker Biomes
fossilSitesClimate2 <- fossilSiteClimate %>%
  filter(!is.na(MAT)) %>%filter(!is.na(MAP))
fossilSitesClimate2 <- whit_it(fossilSitesClimate2, fossilSitesClimate2$MAT, fossilSitesClimate2$MAP) 
filteredModernSites2 <- whit_it(filteredModernSites, filteredModernSites$newMAT, filteredModernSites$newMAP) %>%
  dplyr::select(Site, WhitBiome) %>%
  bind_rows(dplyr::select(fossilSitesClimate2, c(Site, WhitBiome))) %>%
  left_join(logClusters, by = "Site")
allSpecies3 <- allSpecies %>%
  left_join(dplyr::select(filteredModernSites2, c("Site", "WhitBiome")), by = "Site") %>%
  na.omit()

whitFacet <- ggplot(na.omit(filteredModernSites2), aes(y = "", fill = clusterName)) + geom_bar(position = 'fill', color = 'black', width = .1)  +
  facet_wrap(vars(WhitBiome)) + theme_void() + scale_x_reverse() +
  scale_fill_manual(values = palette) + theme(legend.position = "none", strip.text = element_blank())
whitDistFacet <- ggplot() + 
  geom_line(data = allSpecies3, aes(MeanLMA, group = Site, color = clusterName), alpha = 0.5, stat = "density") +
  facet_wrap(vars(WhitBiome)) + scale_x_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  scale_color_manual(values = palette) + theme_bw() + theme(legend.position = "none") + 
  labs(x = expression("Leaf Mass per Area (g/m"^2*")"), y = "Density")
whitFacetLegend <- ggplot(na.omit(filteredModernSites2), aes(y = "", fill = clusterName)) + geom_bar(position = 'fill', color = 'black', width = .1)  +
  facet_wrap(vars(WhitBiome)) + theme_void() + scale_x_reverse() +
  scale_fill_manual("LMA Cluster", values = palette) + 
  theme(legend.position = "right", strip.text = element_blank())
legend <- get_legend(whitFacetLegend)
biomeDistributions <- whitDistFacet + inset_element(whitFacet, 0,0.08,1,1.13)
#Figure6
figure6 <- plot_grid(biomeDistributions, legend, rel_widths = c(1, 0.3))
rm(whitDistFacet, whitFacetLegend, whitFacet, biomeDistributions,legend, allSpecies3, filteredModernSites2)

###########5C: Clusters with deciduousness + species means
decidEval <- filteredModernSites %>%
  dplyr::select(c("Site", "clusterName", "decNum", "evNum", "naNum", "logLMA", "MeanLMA")) %>%
  mutate(clusterName = as.factor(clusterName)) %>%
  mutate(decidP = decNum/(decNum + evNum)*100) %>%
  mutate(lmaBin = cut(logLMA, breaks = 10))

#Tukey test for deciduousness
# fm2 <- aov(decidP ~ clusterName, data = decidEval)
# tuk2 <- glht(fm2, linfct = mcp(clusterName = "Tukey"))
# maptuk <- cld(tuk2)
# maptuk

#Plot deciduousness by cluster
decidCluster <- ggplot(decidEval, mapping = aes(x = clusterName, y = decidP, fill = clusterName)) + geom_boxplot(na.rm = TRUE) + 
  scale_fill_manual(values = palette) + theme_bw() + theme(legend.position = "none") +
  labs(x = "Distribution Cluster", y = "Percent of Deciduous Species at Site") + 
  annotate("rect", xmin = 1.5, xmax = 5.5, ymin = -Inf,ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf,ymax = Inf, alpha = .4) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), text = c("a", "a", "b", "c", "bc", "bc", "c")), aes(x = clusterName, y = -10, label = text)) +
  geom_boxplot(na.rm = TRUE) +
  geom_jitter(width = .1, height = 0, size = .9)

lma4dec <- filteredModernLMA %>% ungroup() %>%
  dplyr::select(Species, SummaryEvDec, MeanLMA, logLMA, clusterName) %>%
  mutate(lmaBin = cut(logLMA, breaks = 10)) %>%
  group_by(lmaBin) %>%
  summarise(nSpecies = n(), nDec = replace_na(length(which(SummaryEvDec=="Deciduous")),0), 
            nEv = replace_na(length(which(SummaryEvDec=="Evergreen")),0), 
            nNA = replace_na(length(which(SummaryEvDec=="NA")),0),
            min = min(MeanLMA),
            max = max(MeanLMA)) %>%
  mutate(chanceofdecid = nDec/(nDec+nEv) * 100)

#Plot decidousness of species per LMA bucket
x <- ggplot(data = lma4dec) + geom_point(aes(y = lmaBin, x = chanceofdecid)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 8.5, alpha = .2) +
  labs(x = "Percent of Species Deciduous", y = expression("Leaf Mass per Area Bin (g/m"^2*")")) + 
  scale_y_discrete(labels = c("(5-8]", "(8-13]", "(13-22]", "(22-36]", "(36-60]", "(60-98]", "(98-162]", "(162-269]", "(269-447]", "(447-741]")) +
  theme_bw() 
y <- ggplot(lma4dec) + geom_bar(aes(y = lmaBin, x = nSpecies), stat = "identity") +
  labs(x = "Species Count") +
  theme_bw() + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank(), 
                     axis.ticks = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank())
decidSpecies <- plot_grid(x,y, align = "h", rel_widths = c(3,1))
#Figure8
figure8 <-  plot_grid(decidSpecies, decidCluster, nrow = 2, labels = c("A", "B"))
rm(decidSpecies, decidCluster, x, y)

#####Part 6: What can we do with the fossil data
##########6A: Strongest correlations with any leaf architectural characters
fossilMorphology2 <- fossilMorphology2 %>%
  dplyr::select(Site, Genus, Species, Source, Lobation, `Margin Type`, `Primary Vein Framework`, `Major Secondary Framework`) %>%
  filter(!is.na(Source)) %>%
  filter(!is.na(`Margin Type`)) %>%
  filter(!is.na(`Primary Vein Framework`)) %>%
  filter(!is.na(`Major Secondary Framework`)) %>%
  filter(!is.na(Lobation))

#Convert toothed untoothed data
marginVect <- data.frame(Site="", Genus="", Species="", Toothed = 0)
for(i in 1:nrow(fossilMorphology2)) {
  marginVect[i,1] <- fossilMorphology2$Site[i]
  marginVect[i,2] <- fossilMorphology2$Genus[i]
  marginVect[i,3] <- fossilMorphology2$Species[i]
  if (fossilMorphology2$`Margin Type`[i]%in% c("Untoothed", "untoothed")) {
    marginVect[i,4] <- 0
  } else if (is.na(fossilMorphology2$`Margin Type`[i])) {
    marginVect[i,4] <- NA
  } else {
    marginVect[i,4] <- 1
  }
  
}

tempFossilLMA <- filteredFossilLMA[c(1:15,62)] %>%
  left_join(fossilMorphology, by = c("Site", "Genus", "Species")) %>%
  filter(!is.na(Source)) %>%
  filter(!is.na(clusterName)) %>%
  filter(ID != "Magnolia dayana") %>%
  filter(ID!= "Trochodendroides HC103") %>%
  left_join(marginVect, by = c("Site", "Genus", "Species")) %>%
  unique()

tempFossilLMA2 <- tempFossilLMA %>%
  dplyr::select(!c("Elliptic", "Obovate", "Ovate", "Oblong", "Linear", "Major Secondary Spacing", "Major Secondary Attachment to Midvein", "Course of percurrent tertiaries", "Angle of percurrent tertiaries", "Intercostal Tertiary Vein Angle Variability", "Quaternary Vein Fabric", "Tooth Spacing", "# Orders of Teeth", "Sinus Shape"))

cols <- c("Medial Symmetry", "Base Symmetry", "Margin Type", "Apex Angle", "Base Angle", "Primary Vein Framework", "Major Secondary Framework", "Intercostal Tertiary Vein Fabric")  
tempFossilLMA2[cols] <- lapply(tempFossilLMA2[cols], factor)  

#Ordinate leaf architectural data
dmat <- gowdis(tempFossilLMA2[18:28], ord = "podani", asym.bin = c(4,11))
ordi <- cmdscale(dmat, eig = TRUE, k = 3)

p1 <- function(){
  par(mar = c(2,3,2,1))
  ordiplot(ordi, choices = c(1,2))
  summary(ordisurf(ordi ~ MeanLMA, tempFossilLMA2, choices = c(1,2), main = "", labcex = 1.5))
  mtext(side = 3, adj = 0, text = "Deviance explained by LMA: 8.15%", cex = 1.5)
}

#Go through every leaf architectural character and compare clusters - Tukey tests are commented out
a <- ggplot(data = subset(tempFossilLMA2, !is.na(`Major Secondary Framework`))) + geom_boxplot(aes(x = (as.character(`Major Secondary Framework`)), y= MeanLMA), fill = "grey") +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  geom_text(data = tibble(`Major Secondary Framework` = c(0:7), text = c("a", "a", "a", "a", "a", "a", "a", "a")), aes(x = (as.character(`Major Secondary Framework`)), y = 350, label = text)) +
  scale_x_discrete(labels = str_wrap(c("Crasp", "Semicrasp", "F.Semicrasp", "Eucampt", "Retic", "Clad", "Brochid", "Mixed"), width = 10)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Major 2° Venation")
#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = `Major Secondary Framework`)
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

b <- ggplot(data = subset(tempFossilLMA2, !is.na(`Medial Symmetry`))) + geom_boxplot(aes(x = as.character(`Medial Symmetry`), y= MeanLMA), fill = "grey")+
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  geom_text(data = tibble(`Medial Symmetry` = c(0:2), text = c("a", "a", "a")), aes(x = (as.character(`Medial Symmetry`)), y = 350, label = text)) +
  scale_x_discrete(labels = str_wrap(c("Sym", "Varies", "Asym"), width = 16)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Medial Symmetry")


#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = `Medial Symmetry`)
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

c <- ggplot(data = subset(tempFossilLMA2, !is.na(`Base Symmetry`))) + geom_boxplot(aes(x = as.character(`Base Symmetry`), y= MeanLMA), fill = "grey")+
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  geom_text(data = tibble(`Base Symmetry` = c(0:2), text = c("a", "a", "a")), aes(x = (as.character(`Base Symmetry`)), y = 350, label = text)) +
  scale_x_discrete(labels = str_wrap(c("Sym", "Varies", "Asym"), width = 16)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Base Symmetry")

#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = `Base Symmetry`)
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

d <- ggplot(data = subset(tempFossilLMA2, !is.na(`Lobation`))) + geom_boxplot(aes(x = as.character(`Lobation`), y= MeanLMA), fill = "grey")+
  geom_text(data = tibble(`Lobation` = c(0:1), text = c("a", "a")), aes(x = (as.character(`Lobation`)), y = 350, label = text)) +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  scale_x_discrete(labels = str_wrap(c("Unlobed", "Lobed"), width = 16)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Lobation")

#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = Lobation) %>% mutate(second = as.factor(second))
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

e <- ggplot(data = subset(tempFossilLMA2, !is.na(`Margin Type`))) + geom_boxplot(aes(x = as.character(`Margin Type`), y= MeanLMA), fill = "grey")+
  geom_text(data = tibble(`Margin Type` = c(0:2), text = c("b", "ab", "a")), aes(x = (as.character(`Margin Type`)), y = 350, label = text)) +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  scale_x_discrete(labels = str_wrap(c("Dentate", "Serrate", "Crenate"), width = 16)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Toothed Margin Type")

#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = `Margin Type`)
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

f <- ggplot(data = subset(tempFossilLMA2, !is.na(`Apex Angle`))) + geom_boxplot(aes(x = as.character(`Apex Angle`), y= MeanLMA), fill = "grey")+
  geom_text(data = tibble(`Apex Angle` = c(0:3), text = c("a", "a", "a", "a")), aes(x = (as.character(`Apex Angle`)), y = 350, label = text)) +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  scale_x_discrete(labels = str_wrap(c("Acute", "Acute to Obtuse", "Obtuse", "Reflex"), width = 10)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Apex Angle")

#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = `Apex Angle`)
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

g <- ggplot(data = subset(tempFossilLMA2, !is.na(`Base Angle`))) + geom_boxplot(aes(x = as.character(`Base Angle`), y= MeanLMA), fill = "grey")+
  geom_text(data = tibble(`Base Angle` = c(0:4), text = c("b", "ab", "a", "ab", "ab")), aes(x = (as.character(`Base Angle`)), y = 350, label = text)) +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  scale_x_discrete(labels = str_wrap(c("Acute", "Acute to Obtuse", "Obtuse", "Obtuse to Reflex" , "Reflex"), width = 10)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Base Angle")

#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = `Base Angle`)
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

#wilcox.test(filter(tempFossilLMA2,`Base Angle`==0)$logLMA, filter(tempFossilLMA2,`Base Angle`==3)$logLMA)
#Acute (0) and Obtuse (3) have sig dif LMA

h <- ggplot(data = subset(tempFossilLMA2, !is.na(`Primary Vein Framework`))) + geom_boxplot(aes(x = as.character(`Primary Vein Framework`), y= MeanLMA), fill = "grey")+
  geom_text(data = tibble(`Primary Vein Framework` = c(0:4), text = c("a", "a", "a", "a", "a")), aes(x = (as.character(`Primary Vein Framework`)), y = 350, label = text)) +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  scale_x_discrete(labels = str_wrap(c("Pinnate", "Actin", "Palinactin", "Basal Acro" , "Suprabasal Acro"), width = 10)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "1° Vein Framework")

#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = `Primary Vein Framework`)
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

i <- ggplot(data = subset(tempFossilLMA2, !is.na(`Intercostal Tertiary Vein Fabric`))) + geom_boxplot(aes(x = as.character(`Intercostal Tertiary Vein Fabric`), y= MeanLMA), fill = "grey")+
  geom_text(data = tibble(`Intercostal Tertiary Vein Fabric` = c(0:2), text = c("a", "a", "a")), aes(x = (as.character(`Intercostal Tertiary Vein Fabric`)), y = 350, label = text)) +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  scale_x_discrete(labels = str_wrap(c("Percurrent", "Reticulate", "Ramified"), width = 10)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Intercostal 3° Vein Fabric")

#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = `Intercostal Tertiary Vein Fabric`)
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

j <- ggplot(data = subset(tempFossilLMA2,!is.na(`Toothed`))) + geom_boxplot(aes(x = as.character(`Toothed`), y= MeanLMA), fill = "grey")+
  geom_text(data = tibble(`Toothed` = c(0:1), text = c("a", "a")), aes(x = (as.character(`Toothed`)), y = 350, label = text)) +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  scale_x_discrete(labels = str_wrap(c("Untoothed", "Toothed"), width = 10)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(),axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Tooth Presence")

#tempe <- tempFossilLMA2 %>% dplyr::select(logLMA, second = Toothed) %>% mutate(second = as.factor(second))
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))

tempFossilLMA3 <- tempFossilLMA2 %>%
  mutate(binRatio = cut(`Laminar L:W Ratio`, c(0, .9, 1.51, 2.51, 3.51, 4.51, 10))) %>%
  subset(!is.na(binRatio))

k <- ggplot(data = (subset(tempFossilLMA3, !is.na(binRatio)))) +  geom_boxplot(aes(x = as.character(binRatio), y = MeanLMA), fill = "grey") +
  geom_text(data = tibble(`binRatio` = levels(tempFossilLMA3$binRatio), text = c("a", "a", "a", "a", "a", "a")), aes(x = (as.character(binRatio)), y = 350, label = text)) +
  #geom_text(data = tibble(`binRatio` = levels(tempFossilLMA3$binRatio), text = as.character(summary(tempFossilLMA3$binRatio, n = n()))), aes(x = (as.character(`binRatio`)), y = 300, label = text)) +
  scale_y_continuous(trans ='log10', breaks = c(10,20,50,100,200,500)) +
  scale_x_discrete(labels = str_wrap(c("<1", "1", "2", "3", "4", "5+"))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(),axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Length:Width Ratio")

#tempe <- tempFossilLMA3 %>% dplyr::select(logLMA, second = binRatio) %>% mutate(second = as.factor(second))
#fm1 <- aov(logLMA ~ second, data = tempe)
#tuk <- glht(fm1, linfct = mcp(second = "Tukey"))
#maptuk <- cld(tuk)
#plot(cld(tuk))


#AppendixS12
appendixS12 <- plot_grid(f,g,b,c,j,e,d,k,h,a,i,NULL, align = "h", axis = "l")
rm(f,g,b,c,j,e,d,k,h,a,i,marginVect)

##########6B: Evaluate TCTs within sites and clusters
#Convert MLA leaf architecutre into TCT assignments
tctList <- vector()
for (i in 1:nrow(fossilMorphology2)) {
  if(fossilMorphology2$Lobation[i] %in% c("Unlobed", "unlobed")) {
    if(fossilMorphology2$`Margin Type`[i] %in% c("Untoothed", "untoothed")) {
      if(fossilMorphology2$`Primary Vein Framework`[i] %in% c("pinnate", "Pinnate")) {
        if(fossilMorphology2$`Major Secondary Framework`[i] %in% c("Brochidodromous", "Mixed", "Semicraspedodromous", "Festooned semicraspedodromous", "Eucamptodromous becoming brochidodromous distally","Simple brochidodromous", "Simple Brochidodromous", "Festooned brochidodromous")) {
          tctList[i] <- "A"
        }
        else {
          tctList[i] <- "B"
        }
      }
      else {
        if(fossilMorphology2$`Major Secondary Framework`[i] %in% c("Brochidodromous", "Mixed", "Semicraspedodromous", "Festooned semicraspedodromous", "Eucamptodromous becoming brochidodromous distally","Simple brochidodromous", "Simple Brochidodromous", "Festooned brochidodromous")) {
          tctList[i] <- "C"
        }
        else {
          tctList[i] <- "D"
        }
      }
    }
    else {
      if(fossilMorphology2$`Primary Vein Framework`[i] %in% c("pinnate", "Pinnate")) {
        if(fossilMorphology2$`Major Secondary Framework`[i] %in% c("Brochidodromous", "Mixed", "Semicraspedodromous", "Festooned semicraspedodromous", "Eucamptodromous becoming brochidodromous distally","Simple brochidodromous", "Simple Brochidodromous", "Festooned brochidodromous")) {
          tctList[i] <- "E" 
        }
        else {
          tctList[i] <- "F"
        }
      }
      else {
        if(fossilMorphology2$`Major Secondary Framework`[i] %in% c("Brochidodromous", "Mixed", "Semicraspedodromous", "Festooned semicraspedodromous", "Eucamptodromous becoming brochidodromous distally","Simple brochidodromous", "Simple Brochidodromous", "Festooned brochidodromous")) {
          tctList[i] <- "G"
        }
        else {
          tctList[i] <- "H"
        }
      }
    }
  }
  else{
    if(fossilMorphology2$`Margin Type`[i] %in% c("Untoothed", "untoothed")) {
      if(fossilMorphology2$`Primary Vein Framework`[i] %in% c("pinnate", "Pinnate")) {
        if(fossilMorphology2$`Major Secondary Framework`[i] %in% c("Brochidodromous", "Mixed", "Semicraspedodromous", "Festooned semicraspedodromous", "Eucamptodromous becoming brochidodromous distally","Simple brochidodromous", "Simple Brochidodromous", "Festooned brochidodromous")) {
          tctList[i] <- "I"
        }
        else {
          tctList[i] <- "J"
        }
      }
      else {
        if(fossilMorphology2$`Major Secondary Framework`[i] %in% c("Brochidodromous", "Mixed", "Semicraspedodromous", "Festooned semicraspedodromous", "Eucamptodromous becoming brochidodromous distally","Simple brochidodromous", "Simple Brochidodromous", "Festooned brochidodromous")) {
          tctList[i] <- "K"
        }
        else {
          tctList[i] <- "L"
        }
      }
    }
    else {
      if(fossilMorphology2$`Primary Vein Framework`[i] %in% c("pinnate", "Pinnate")) {
        if(fossilMorphology2$`Major Secondary Framework`[i] %in% c("Brochidodromous", "Mixed", "Semicraspedodromous", "Festooned semicraspedodromous", "Eucamptodromous becoming brochidodromous distally","Simple brochidodromous", "Simple Brochidodromous", "Festooned brochidodromous")) {
          tctList[i] <- "M" 
        }
        else {
          tctList[i] <- "N"
        }
      }
      else {
        if(fossilMorphology2$`Major Secondary Framework`[i] %in% c("Brochidodromous", "Mixed", "Semicraspedodromous", "Festooned semicraspedodromous", "Eucamptodromous becoming brochidodromous distally","Simple brochidodromous", "Simple Brochidodromous", "Festooned brochidodromous")) {
          tctList[i] <- "O"
        }
        else {
          tctList[i] <- "P"
        }
      }
    }
  }
}

neoMorphoTCT <- fossilMorphology2 %>%
  bind_cols(TCT =tctList) %>%
  dplyr::select(Group = Site, Genus, Species, TCT)

filteredFossilTCT <- filteredFossilLMA %>%
  mutate(Group = Site)
filteredFossilTCT$Group[143:389] <-  "Denver Basin"

neoMorphoTCT2 <- filteredFossilTCT %>%
  left_join(neoMorphoTCT, by = c("Group", "Genus", "Species")) %>%
  ungroup() %>%
  unique()

finalTCTDataset <- neoMorphoTCT2 %>% dplyr::select(logLMA, TCT) %>% mutate(TCT = as.factor(TCT)) %>%
  filter(TCT %in% c("A", "B", "C", "D", "E", "F", "G", "H", "K", "N", "O", "P"))
#fm1 <- aov(logLMA ~ TCT, data = finalTCTDataset)
#tuk <- glht(fm1, linfct = mcp(TCT = "Tukey"))
#maptuk <- cld(tuk)
#par(mar = c(5,4,4,2))

#AppendixS13
appendixS13 <- ggplot(data = finalTCTDataset) + geom_boxplot(aes(x = as.character(TCT), y= 10^logLMA), fill = "grey")+
  geom_point(aes(x = as.character(TCT), y = 10^logLMA)) +
  geom_text(data = tibble(`TCT` = c("A", "B", "C", "D", "E", "F", "G", "H", "K", "N", "O", "P"), 
                          text = c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a")), 
            aes(x = (as.character(TCT)), y = 350, label = text)) +
  scale_y_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) +
  
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x = element_blank(),axis.title.y = element_text(), plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "TCT", y = "LMA")

TCTSite <- neoMorphoTCT2 %>%
  dplyr::select(Site, TCT, clusterName, n) %>%
  filter(!is.na(clusterName)) %>%
  mutate(count = 1) %>%
  pivot_wider(names_from = TCT, values_from = count, values_fn = sum)

TCTSite2 <- TCTSite %>%
  filter(`NA`/n < 0.5) %>%
  left_join(dplyr::select(filteredFossilSites, c(Site, MeanLMA, logLMA)), by = "Site")

rownames(TCTSite2) <- TCTSite2$Site

TCTSite3 <- TCTSite2 %>%
  dplyr::select(!c(I,J,M, `NA`, `clusterName`, n, MeanLMA, logLMA)) %>%
  dplyr::select(!Site) %>%
  mutate_all(~replace_na(.,0))

rownames(TCTSite3) <- TCTSite2$Site

#NMDS with TCT data at each site, then compare by clusters
tctNMDS <- metaMDS(TCTSite3, k =3, distance = "bray", trymax =  150)
tct_scores <- as.data.frame(scores(tctNMDS)$sites) %>%
  rownames_to_column(var = "Site") %>%
  left_join(TCTSite2, by = "Site") %>%
  dplyr::select(Site, NMDS1, NMDS2, clusterName)
tct_species <- as.data.frame(scores(tctNMDS)$species)%>%
  rownames_to_column(var= "TCT") %>%
  mutate(Toothed = c(1,1,1,0,0,1,1,1,1,0,0,0,0)) %>%
  mutate(Lobed = c(0,0,0,0,0,1,1,0,1,0,1,0,1)) %>%
  mutate(Pinnate = c(1,1,0,1,0,1,0,0,0,1,0,0,0)) %>%
  mutate(Looped = c(1,0,1,0,1,0,0,0,1,1,0,0,1))

p2 <- function() {
  par(mar = c(2,3,2,1))
  ordiplot(tctNMDS, choices = c(1,2), display = "sites")
  ordispider(tctNMDS, choices =c(1,2), groups = TCTSite2$clusterName, label = FALSE, col = palette[2:5])
  text(tctNMDS$species, rownames(tctNMDS$species), col = "red")
  legend("bottomleft", col = palette[3:5], lty =1, bty ="n", legend = c("Middle Left", "Middle Right", "Low Right"))
}

#Anosim for clusters with TCTs
attach(TCTSite2)
tctNMDS.ano <- anosim(TCTSite3, grouping = clusterName)
#plot(tctNMDS.ano)
#summary(tctNMDS.ano)
detach(TCTSite2)

#Figure9
figure9 <- plot_grid(p1, p2, nrow = 2, align = "hv", labels = c("A", "B"))
rm(p1, p2)

##########6C: Clusters with Depositional Environment
filteredFossilSites$DepEnv <- factor(filteredFossilSites$DepEnv, levels = c("Fluvial", "Coastal", "Lacustrine", "Mixed"))
filteredFossilLMA$DepEnv <- factor(filteredFossilLMA$DepEnv, levels = c("Fluvial", "Coastal", "Lacustrine", "Mixed"))

depBarPlot <- ggplot(filteredFossilSites) + geom_bar(aes(y = "", fill = clusterName), position = 'fill', width = .1, color = 'black') +
  facet_wrap(vars(DepEnv)) + scale_fill_manual(values = palette[2:6]) + scale_x_reverse() +
  theme_void() + theme(legend.position = "none", strip.text = element_blank()) 
depLinePlot <- ggplot(filter(filteredFossilLMA, !is.na(clusterName))) + geom_line(aes(x = MeanLMA, group = Site, color = clusterName), alpha = 0.5, stat = "density") +
  facet_wrap(vars(DepEnv)) + scale_x_continuous(trans = 'log10', breaks = c(10,20,50,100,200,500)) + scale_color_manual(values = palette[2:6]) +
  theme_bw() + theme(legend.position = "none") + labs(x = expression("Leaf Mass per Area (g/m"^2*")"), y = "Density")
depBarPlotLeg <- ggplot(filteredFossilSites) + geom_bar(aes(y = "", fill = clusterName), position = 'fill', width = .1, color = 'black') +
  facet_wrap(vars(DepEnv)) + scale_fill_manual("LMA Cluster", values = palette[2:6]) + scale_x_reverse() +
  theme_void() + theme(legend.position = "right", strip.text = element_blank()) 
deplegend <- get_legend(depBarPlotLeg)
bigDepPlot <- depLinePlot+ inset_element(depBarPlot, 0,0.15,1,1.2)
#Figure10
figure10 <- plot_grid(bigDepPlot, deplegend, rel_widths = c(1, 0.3))
rm(bigDepPlot, deplegend, depBarPlotLeg, depLinePlot, depBarPlot)

#Test whether there are differences in depositional environment between clusters
depEnvdata <- filteredFossilSites %>%
  dplyr::select(c("Site", "DepEnv", "clusterName", "logLMA", "MeanLMA")) %>%
  filter(DepEnv == "Fluvial" | DepEnv == "Lacustrine") %>% 
  group_by(DepEnv) %>%
  summarise(meanLMA = mean(MeanLMA), logLMA = mean(logLMA))
datDot <- matrix(c(9,17,9,9,0,1,4,3), ncol = 4, nrow =2)
rownames(datDot) <- c("Lacustrine", "Fluvial")
colnames(datDot) <- c("Middle Left", "Middle Right", "Low Left", "Low Right")
datDot <- as.table(datDot)
depositionalChiSquare <- chisq.test(datDot)

################7: Clusters with environmental conditions
filteredModernSites11 <- filteredModernSites
filteredModernSites11$clusterName <- as.factor(filteredModernSites$clusterName)
filteredFossilSites11 <- filteredFossilSites
filteredFossilSites11$clusterName <- as.factor(filteredFossilSites$clusterName)

#fm1 <- aov(log10(newMAP) ~ clusterName, data = filteredModernSites11)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

mapCluster <- ggplot(filteredModernSites, aes(x = clusterName, y = newMAP, fill = clusterName)) + 
  geom_boxplot(notch = FALSE, outlier.alpha = 0) + 
  annotate("rect", xmin = 1.5, xmax = 5.5, ymin = -Inf,ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf,ymax = Inf, alpha = .4) +
  geom_boxplot(notch = FALSE, outlier.alpha = 0) +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), text = c("a", "ab", "ac", "d", "c", "acd", "bcd")), aes(x = clusterName, y = -1, label = text)) +
  scale_fill_manual(values = palette) +
  #scale_y_continuous(trans = 'log10') +
  labs(x = "", y = "Mean Annual Precipitation (mm)") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_blank())

#fm1 <- aov(newMAT ~ clusterName, data = filteredModernSites11)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

matCluster <- ggplot(filteredModernSites, aes(x = clusterName, y = newMAT, fill = clusterName)) +
  geom_boxplot(notch = FALSE, outlier.alpha = 0) + 
  annotate("rect", xmin = 1.5, xmax = 5.5, ymin = -Inf,ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf,ymax = Inf, alpha = .4) +
  geom_boxplot(notch = FALSE, outlier.alpha = 0) +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), text = c("a", "a", "b", "c", "bc", "abc", "bc")), aes(x = clusterName, y = -10, label = text)) +
  scale_fill_manual(values = palette) +
  labs(x = "", y = expression("Mean Annual Temperature ("*degree*C*")")) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_blank())

#fm1 <- aov(gridVPD ~ clusterName, data = filteredModernSites11)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

vpdCluster <- ggplot(filteredModernSites, aes(x = clusterName, y = gridVPD, fill = clusterName)) + 
  geom_boxplot(notch = FALSE, outlier.alpha = 0) + 
  annotate("rect", xmin = 1.5, xmax = 5.5, ymin = -Inf,ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf,ymax = Inf, alpha = .4) +
  geom_boxplot(notch = FALSE, outlier.alpha = 0) +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), text = c("a", "ab", "bc", "c", "c", "ac", "c")), aes(x = clusterName, y = -1, label = text)) +
  scale_fill_manual(values = palette) +
  labs(x = "", y = "Vapor Pressure Density (Pa)") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

#fm1 <- aov(gridPET ~ clusterName, data = filteredModernSites11)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

petCluster <- ggplot(filteredModernSites, aes(x = clusterName, y = gridPET, fill = clusterName)) + 
  geom_boxplot(notch = FALSE, outlier.alpha = 0) + 
  annotate("rect", xmin = 1.5, xmax = 5.5, ymin = -Inf,ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf,ymax = Inf, alpha = .4) +
  geom_boxplot(notch = FALSE, outlier.alpha = 0) +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), text = c("a", "a", "b", "b", "b", "ab", "b")), aes(x = clusterName, y = -1, label = text)) +
  scale_fill_manual(values = palette) +
  labs(x = "", y = expression("Potential Evapotranpsiration (kg/m"^2*")")) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle =90, hjust = 1, vjust = .5))

#fm1 <- aov(gridRSDS ~ clusterName, data = filteredModernSites11)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

rsdsCluster <- ggplot(filteredModernSites, aes(x = clusterName, y = gridRSDS, fill = clusterName)) + 
  geom_boxplot(notch = FALSE, outlier.alpha = 0) + 
  annotate("rect", xmin = 1.5, xmax = 5.5, ymin = -Inf,ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf,ymax = Inf, alpha = .4) +
  geom_boxplot(notch = FALSE, outlier.alpha = 0) +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), text = c("a", "a", "ab", "bc", "bc", "ac", "c")), aes(x = clusterName, y = -1, label = text)) +
  scale_fill_manual(values = palette) +
  labs(x = "", y = expression("Surface Downwelling Shortwave Radiation (MJ/m"^2*"d)")) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

#fm1 <- aov(TempSeasonality ~ clusterName, data = filteredModernSites11)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

TempSeasonCluster <- ggplot(filteredModernSites, aes(x = clusterName, y = TempSeasonality, fill = clusterName)) + 
  geom_boxplot(notch = FALSE, outlier.alpha = 0) + 
  annotate("rect", xmin = 1.5, xmax = 5.5, ymin = -Inf,ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf,ymax = Inf, alpha = .4) +
  geom_boxplot(notch = FALSE, outlier.alpha = 0) +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), text = c("a", "b", "c", "d", "cd", "bcd", "cd")), aes(x = clusterName, y = -1, label = text)) +
  scale_fill_manual(values = palette) +
  labs(x = "", y = expression("Temperature Seasonality (SD * 100)")) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

#fm1 <- aov(co2 ~ clusterName, data = filteredFossilSites)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

co2Cluster <- ggplot(subset(filteredFossilSites, !is.na(clusterName)), aes(x = clusterName, y = co2, fill = clusterName)) + 
  geom_boxplot(notch = FALSE, outlier.alpha = 0) + 
  annotate("rect", xmin = .5, xmax = 4.5, ymin = -Inf,ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 1.5, xmax = 3.5, ymin = -Inf,ymax = Inf, alpha = .4) +
  annotate("rect", xmin = -.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha = 0) +
  geom_boxplot(notch = FALSE, outlier.alpha = 0) +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right"), text = c("a", "a", "a", "a")), aes(x = clusterName, y = -1, label = text)) +
  scale_fill_manual(values = palette[2:5]) +
  labs(x = "", y = expression("Atmospheric CO"[2]*" (ppm)")) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_blank())

#Figure5
figure5 <- matCluster + mapCluster + co2Cluster + rsdsCluster + TempSeasonCluster + vpdCluster + plot_annotation(tag_levels = 'A')
rm(matCluster, mapCluster, co2Cluster, rsdsCluster, TempSeasonCluster, vpdCluster)

filteredModernSites3 <- filteredModernSites %>%
  dplyr::select(c(2, 49, 12:21, 23:29, 32:35, 37:43, 46:47)) %>%
  na.omit() %>%
  mutate(newMAP = log10(newMAP), PrecipWettestMonth = log10(PrecipWettestMonth), PrecipDriestMonth = log10(PrecipDriestMonth+1), PrecipWettestQuarter = log10(PrecipWettestQuarter), PrecipDriestQuarter = log10(PrecipDriestQuarter), PrecipWarmestQuarter = log10(PrecipWarmestQuarter), PrecipColdestQuarter = log10(PrecipColdestQuarter))

#AppendixS11
pca <- rda(filteredModernSites3[3:32], scale = TRUE)
par(cex = 1)
biplot(pca, display = c("sites", "species"), type = c("text", "points"))
legend("topright", col = palette, lty =1, legend = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"))
ordiellipse(pca, group = filteredModernSites3$clusterName, col = palette, label = FALSE, kind = "sd")
appendixS11 <- recordPlot()

################7: Evaluate variance of clusters
variance <- filteredModernLMA %>%
  group_by(Database, Site, clusterName) %>%
  summarize(variance = (var(logLMA)))
variance$clusterName <- as.factor(variance$clusterName)
#fm1 <- aov(variance ~ clusterName, data = variance)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

clusterVarPlot <- ggplot(variance, aes(x = clusterName, y = variance, fill = clusterName)) +
  geom_boxplot() +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), 
                              text = c("ab", "a", "b", "b", "a", "b", "b")), aes(x = clusterName, y = 0, label = text)) +
  scale_fill_manual(values = palette) +
  labs(x = "", y = "Site Variance of species log10(Leaf Mass per Area)") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


#filteredModernSites$clusterName <- as.factor(filteredModernSites$clusterName)
#fm1 <- aov(logLMA ~ clusterName, data = filteredModernSites)
#tuk <- glht(fm1, linfct = mcp(clusterName = "Tukey"))
#maptuk <- cld(tuk)
#maptuk

clusterMeanPlot <- ggplot(filteredModernSites, aes(x = clusterName, y = logLMA, fill = clusterName)) +
  geom_boxplot() +
  geom_jitter(width = .1, height = 0, size = .9) +
  geom_text(data = data.frame(clusterName = c("1. Middle Leftmost", "2. Low Left", "3. Middle Left", "4. Middle Right", "5. Low Right", "6. Upper Right", "7. Middle Rightmost"), 
                              text = c("a", "b", "c", "d", "e", "f", "f")), aes(x = clusterName, y = 1, label = text)) +
  scale_fill_manual(values = palette) +
  labs(x = "", y = "Site Mean of species log10(Leaf Mass per Area)") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

#AppendixS10
appendixS10 <- clusterMeanPlot + clusterVarPlot + plot_annotation(tag_levels = 'A')
rm(clusterMeanPlot, clusterVarPlot)

#########Extra - print tables for appendices
#Print some tables
appendixS3 <- filteredFossilLMA %>%
  dplyr::select("Site", "ID", "MeanLMA", "samples") %>%
  left_join(data.frame(Site = neoMorphoTCT2$Site, ID = neoMorphoTCT2$ID, TCT = neoMorphoTCT2$TCT), by = c("Site", "ID")) %>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

appendixS4 <- filteredFossilSites %>%
  left_join(data.frame(Site = fossilSitesClimate2$Site, WhitBiome= fossilSitesClimate2$WhitBiome), by = "Site") %>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

appendixS5 <- filteredModernLMA %>%
  dplyr::select("Dataset" = Database, Site, Species, "Leaf Habit" = SummaryEvDec, "Mean LMA" = MeanLMA, n) %>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

appendixS6 <- whit_it(filteredModernSites, filteredModernSites$newMAT, filteredModernSites$newMAP) %>%
  left_join(dplyr::select(variance, c(Database, Site, variance)), by = c("Database","Site")) %>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 4))) %>%
  mutate(across(c("Lat", "Lon", "MeanLMA", ), \(x) round(x, digits = 2)))

appendixS1 <- tempFossilLMA2 %>%
  dplyr::select(-Family, -Genus, -Species, -samples, -logLMA, -n, -Age, -DepEnv, -MAT, -MAP, -TempBin, -PrecipBin, -clusterName) %>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))
