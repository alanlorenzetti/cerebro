# alorenzetti 20210120

# description ####
# this script will take
# results of cerebro
# and parse to generate
# results

# loading libraries ####
if(!require(pacman)){install.packages("pacman")}
library(pacman)

# required packages
packs = c("tidyverse",
          "ggpubr")

# loading
p_load(char = packs)

# setting ggtheme
theme_set(theme_bw())

# reading features table ####
# features.txt described in README
featureTable = read_tsv("resultsR-area100-circ08-eccen05/features.txt")

# learning thresholds from feature table
# and classifying colonies
intensityThr = quantile(featureTable$b.mean, probs = 0.99)
intensitySdThr = quantile(featureTable$b.sd, probs = 0.01)
mutantIndex = featureTable$b.mean >= intensityThr & featureTable$b.sd <= intensitySdThr
featureTable$class = NA
featureTable$class[mutantIndex] = "Mutant"
featureTable$class[!mutantIndex] = "Normal"
#clusters = featureTable %>%
#  select(b.mean, b.sd) %>%
#  kmeans(., centers = 2, iter.max = 100)
#featureTable$kmeans = clusters$cluster

# plotting scatter with decision boundaries
decisionplot = featureTable %>% 
  ggplot(aes(x=b.mean, y=b.sd, colour=class)) +
  geom_point(alpha=0.2) +
  scale_colour_manual(values = c("Mutant" = "#E15759", "Normal" = "#59A14F"),
                      labels = c("Mutant" = "Mutante", "Normal" = "Normal"),
                      guide = guide_legend(title = "Classe")) +
  geom_vline(xintercept = intensityThr) +
  geom_hline(yintercept = intensitySdThr) +
  xlab("Intensidade Média da Cor Vermelha") +
  ylab("Desvio Padrão da Intensidade da Cor Vermelha")

# getting general knowledge about frequencies
general = featureTable %>% 
  mutate(plateID = str_replace(plate, "^.*/(.*).tif$", "\\1")) %>% 
  dplyr::select(plateID,
                class) %>% 
  group_by(plateID) %>% 
  summarise(counts = n(),
            mutcounts = sum(class == "Mutant")) %>% 
  ungroup() %>% 
  mutate(strain = str_replace(plateID, "^(.*)_.*_.*$", "\\1"),
         replicate = str_replace(plateID, "^.*_(.*_.*)$", "\\1")) %>% 
  group_by(strain) %>% 
  summarise(nPlates = length(strain),
            counts = sum(counts),
            mutcounts = sum(mutcounts),
            frequency = sum(mutcounts)/sum(counts)) %>% 
  mutate(strain = case_when(strain == "dLSm" ~ "dUra3dLSm",
                            strain == "dUra3" ~ "dUra3",
                            strain == "NRC1" ~ "NRC-1",
                            TRUE ~ as.character(strain)))

# plotting general insights
generalplot = general %>%
  ggplot(aes(x=strain, y=frequency,
             label=paste0("p: ", nPlates, "\nc: ", counts))) +
  geom_bar(stat = "identity", fill="white", colour="black") +
  geom_text(vjust=-0.3, size=4) +
  scale_x_discrete(limits = c("NRC-1", "dUra3", "dUra3dLSm")) +
  scale_y_continuous(limits = c(0, 0.025)) +
  labs(x = "Linhagem", y = "Frequência de Mutantes Translúcidos")

# arranging plots
arranged = ggarrange(plotlist = list(decisionplot, generalplot),
                     labels = "AUTO",
                     widths = c(1.5,1))

# saving
ggsave(filename = "panelEspontaneousMutants.png",
       plot = arranged,
       width = 8,
       height = 5,
       units = "in")

# performing Fisher Exact test ####
# to test if strains have different frequencies
# dura3 vs nrc1
dura3Mutants = general[2,4]
nrc1Mutants = general[3,4]

dura3Normal = general[2,3] - general[2,4]
nrc1Normal = general[3,3] - general[3,4]

M = matrix(c(dura3Normal, nrc1Normal, dura3Mutants, nrc1Mutants) %>% unlist() %>% unname(),
           nrow=2,
           dimnames = list(Genotype=c("dUra3", "NRC1"),
                           Phenotype=c("Normal", "Mutant")))

FisherTest = fisher.test(M, alternative = "two.sided")
pvalue = FisherTest$p.value

# dLSm vs nrc1
dLSmMutants = general[1,4]
nrc1Mutants = general[3,4]

dLSmNormal = general[1,3] - general[1,4]
nrc1Normal = general[3,3] - general[3,4]

M = matrix(c(dLSmNormal, nrc1Normal, dLSmMutants, nrc1Mutants) %>% unlist() %>% unname(),
           nrow=2,
           dimnames = list(Genotype=c("dLSm", "NRC1"),
                           Phenotype=c("Normal", "Mutant")))

FisherTest = fisher.test(M, alternative = "two.sided")
pvalue = FisherTest$p.value

# dLSm vs dura3
dLSmMutants = general[1,4]
dura3Mutants = general[2,4]

dLSmNormal = general[1,3] - general[1,4]
dura3Normal = general[2,3] - general[2,4]

M = matrix(c(dLSmNormal, dura3Normal, dLSmMutants, dura3Mutants) %>% unlist() %>% unname(),
           nrow=2,
           dimnames = list(Genotype=c("dLSm", "dUra3"),
                           Phenotype=c("Normal", "Mutant")))

FisherTest = fisher.test(M, alternative = "two.sided")
pvalue = FisherTest$p.value
