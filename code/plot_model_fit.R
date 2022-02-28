library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(lemon)
library(cowplot)
library(ggpubr)
source("./functions_auxiliary.R")


regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1
dataAlpha <- 0.5


############################
###
### Plot simulated data
###
############################
simData <- read.csv("../data/raw_data/simulated_data.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., infections1=Prevalence1*Population/100,
                prevChange=Prevalence2-Prevalence1,
                infections2=prevChange*Population/100,
                deathChange=Deaths2-Deaths1,
                meanAge=mid_bin_age(Age),
                Location=factor(Location)) %>%
  dplyr::filter(., Prevalence1>0)

simPlot <- ggplot(simData, aes(x=meanAge, y=Deaths1/infections1*100, group=Location)) +
  geom_point(size=locPointSize, alpha=dataAlpha, color="#1B9E77") +
  geom_line(size=locLineSize, alpha=dataAlpha, color="#1B9E77") +
  geom_point(aes(x=meanAge, deathChange/infections2*100), color="#E7298A",
             size=locPointSize, alpha=dataAlpha) +
  geom_line(aes(x=meanAge, deathChange/infections2*100), color="#E7298A",
             size=locLineSize, alpha=dataAlpha) +
  theme_bw() +
  scale_y_continuous(trans='log10', labels=scaleFun,
                     breaks=10^c(-3, -2, -1, 0, 1, 2)) +
  theme(strip.background=element_rect(fill="white", color="white"),
        strip.text=element_text(face="bold"),
        legend.position="top",
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour = "black"),
        axis.line.x=element_line(size=0.5, linetype="solid"),
        axis.line.y=element_line(size=0.5, linetype="solid"),
        legend.title=element_blank()) +
  xlab("Age") +
  ylab("% Deaths")

ggsave("../data/figures/simulated_data.png", simPlot, units="cm",
       height=13, width=15)


posteriorTraces <- read.csv("../data/analysis_results/simulated_data_fit.csv") %>%
  as_tibble(.)

serologyTrace <- dplyr::mutate(posteriorTraces, .chain=factor(.chain))  %>%
  dplyr::filter(., !(.variable %in% c("locationIntercept", "locationChange"))) %>%
  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
  geom_line() +
  facet_wrap(.~.variable, scales="free") +
  theme_bw()

ggsave("../data/figures/simulated_data_traces.png", serologyTrace, units="cm",
       height=15, width=22)

parameters <- c("ageSlope", "intercept", "interceptChange", "interceptSigma",
                "changeSigma")
parameterVals <- c(0.13*26.13, -12.9, -0.7, 1.11, 0.1) 
paramPlots <- list()
for (p in c(1:length(parameters))) {
  paramPlots[[p]] <- dplyr::filter(posteriorTraces, .variable==parameters[p]) %>%
    ggplot(., aes(x=.value)) +
    geom_density(color="lightblue", fill="lightblue") +
    geom_vline(xintercept=parameterVals[p]) +
    ggtitle(parameters[p])
}

allParams <- ggpubr::ggarrange(plotlist=paramPlots)

ggsave("../data/figures/simulated_data_fit_estimates.png", allParams,
       units="cm", height=10, width=22)


