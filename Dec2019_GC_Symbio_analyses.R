#Symbiodiniaceae community diversity analyses for Genetic Contraints project
#edited on 12.16.19, Lauren Howe-Kerr

# load libraries
library('phyloseq'); packageVersion('phyloseq') #version 1.26.1
library('vegan'); packageVersion('vegan') #version 2.5.6
library('dplyr'); packageVersion('dplyr') #version 0.8.3
library('ggplot2'); packageVersion('ggplot2') #version 3.2.1
library(DESeq2); packageVersion('DESeq2') #version 1.22.2

#testing diversity analyses with phyloseq object made from LULU ran with 84% and 70% sequence similarity and co-occurrence parameters, respectively
GC <- readRDS("/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/R2/LULUDec19GCiterations/LULUps_GC_Dec2019.rds")

###########################################################
##############Alpha diversity##############################
###########################################################
#now using phyloseq object made from LULU ran with 84% and 95% sequence similarity and co-occurrence parameters, respectively
GC <- readRDS("/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/Analysis_resubmission_GC_Aug2019/LULUps_GC_Aug2019.rds")
hist(sample_sums(GC))
GC <- prune_samples(sample_sums(GC)>=10000, GC) #got rid of 2 samples (G4C6, G4H23)
GC <- subset_taxa(GC, (Class!="NA")) #got rid of one weird ASV that doesnt match to any symbiodiniaceae type on NCBI

ps <- GC
qd <- subset_samples(ps, Treatment2 =="control") #just want to look at communities under control conditions


sample_data(qd)$Genotype <- factor(sample_data(qd)$Genotype, levels= c("G4", "G12", "G20", "G38", "G2", "G27", "G31", "G34"))

levels(sample_data(qd)$Genotype) <-  c("Best-4", "Best-12", "Best-20", "Best-38", "Worst-2", "Worst-27", "Worst-31", "Worst-34")

alphaplot <- plot_richness(qd, x="Genotype", measures=c("Shannon", "Simpson"), color="Status") + 
  labs(x="Genet") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.background = element_rect(fill="#FFFFFF"))
alphaplot

alphaplot_color <- alphaplot + scale_color_manual(values=c("#4169E1", "firebrick"))
ggsave("/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/R2/controlalphadiv_dec19.jpeg", width=6, height=4)

#stats
erich <- estimate_richness(qd, measures=c("Shannon", "Simpson"))
#make dataframe
SampleID <- row.names(erich) 
# makes a new dataframe with the specified columnes
alpha <- data.frame(SampleID, erich) 
s <- data.frame(sample_data(qd)) 
# Change first column title to SampleID to match distances dataframe
colnames(s)[1] <- "SampleID" 
# merges metadata df and distances df
alphaGC <- merge(alpha, s, by = "SampleID") 

#shannon diversity- best 0.8770933, worst 1.2069904
#simpson diversity- best 0.4726877, worst 0.6514613
tapply(alphaGC$Shannon, alphaGC$Status, mean)
tapply(alphaGC$Simpson, alphaGC$Status, mean)

tapply(alphaGC$Shannon, alphaGC$Genotype, mean)
tapply(alphaGC$Simpson, alphaGC$Genotype, mean)

#test for normal distribution
shapiro.test(alphaGC$Shannon) #not normal 
shapiro.test(alphaGC$Simpson) #not normal

wilcox_stat <- wilcox.test(Simpson ~ Status, data=alphaGC)
wilcox_stat

wilcox_stat <- wilcox.test(Shannon ~ Status, data=alphaGC)
wilcox_stat

###########################################################
############Transforming data##############################
###########################################################

#first transform counts data using DESeq2
qduntransformed <- qd
diagdds = phyloseq_to_deseq2(qd, ~ Genotype)

#rlog transform
rlogCOUNTS<-rlog(diagdds,blind=TRUE)
head(assay(rlogCOUNTS))

dat=as.data.frame(assay(rlogCOUNTS)) 

dat[dat < 0.0] <- 0.0 #remove negative values before calculating bray curtis distances
otu_table(qd) <- otu_table(dat, taxa_are_rows = TRUE) #put this transformed data back into phyloseq object as OTU table

write.csv(otu_table(qd), file = "/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/R2/deseqcountsGCDec19.csv")

#qd <- prune_taxa(taxa_sums(qd) > 0, qd) 

#proceed with beta diversity analyses below with this new transformed qd phyloseq object

###########################################################
############NMDS & beta div analysis#######################
###########################################################
GC_bc <- phyloseq::distance(qd, method = "bray")

ord_bc <- ordinate(qd, "NMDS", distance = GC_bc)

statusplot <- plot_ordination(qd, ord_bc, color = "Status")

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls <- plyr::ddply(statusplot$data, "Status", find_hull)

nmdsdata <- plot_ordination(qd, ord_bc)$data
statusplotborder <- ggplot(nmdsdata, aes(x = NMDS1, y = NMDS2, col = factor(Status), fill = factor(Status)))+
  geom_polygon(data = hulls, alpha = 0.2) + geom_point(aes(col = Status, fill = Status)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))

betaplot <- statusplotborder + scale_color_manual(values=c("#4169E1", "firebrick")) + scale_fill_manual(values=c("#4169E1", "firebrick"))

ggarrange(alphaplot, statusplotborder, 
          labels = c("A", "B"),
          ncol = 2,
          legend = "right", common.legend = TRUE)

ggsave("/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/R2/host_div_combo_fig3.jpeg", width=9, height=4)
ggsave("/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/Re2/host_div_combo_fig3_sadcolors.jpeg", width=9, height=4)

# Use bray curtis distance matrix (GC_bc) from above

# make a data frame from the sample_data
qd_sampledf <- data.frame(sample_data(qd))

# Adonis test
adonis(GC_bc ~ Status, data = qd_sampledf)

# Homogeneity of dispersion test
beta <- betadisper(GC_bc, qd_sampledf$Status)
permutest(beta)

###########################################################
####taxonomic assignment/rel abundance bargraph############
###########################################################
library(colorspace)
NewSampleNames.ps <- GC

allotus <- names(sort(taxa_sums(NewSampleNames.ps), decreasing=TRUE))
ps.all <- transform_sample_counts(NewSampleNames.ps, function(OTU) OTU/sum(OTU)) #relative abundance by sample

ps.all.over.1 <- filter_taxa(ps.all, function(OTU) mean(OTU) > .0001, TRUE) #removing low abundance samples
allotus.over.1 <- names(sort(taxa_sums(ps.all.over.1), decreasing=TRUE))
ps.all.over.1 <- prune_taxa(allotus.over.1, ps.all.over.1)

as.data.frame(otu_table(ps.all))

NewSampleNames.ps <- ps.all.over.1
sample_data(NewSampleNames.ps)$Sample <- c("All1","All2",	"All3",	"Bac1",	"Bac2",	"Bac3",	"Cont1",	"Cont2",	"Cont3",	"Heat1",	"Heat2", "pCO2-1",	"pCO2-2",	"pCO2-3",	"All1",	"All2",	"All3",	"Bac1",	"Bac2",	"Bac3",	"Cont1",	"Cont2",	"Cont3",	"Heat1",	"Heat2",	"Heat3",	"pCO2-1",	"pCO2-2",	"pCO2-3",	"All1",	"All2",	"Bac1",	"Bac2",	"Bac3",	"Cont1",	"Cont2",	"Cont3",	"Heat1",	"pCO2-1",	"Bac1",	"Bac2",	"Cont1",	"Cont2",	"Cont3",	"pCO2-1",	"All1",	"All2",	"Bac1",	"Bac2",	"Cont1",	"Cont2",	"Cont3",	"Heat1",	"Heat2",	"pCO2-1",	"pCO2-2",	"pCO2-3",	"All1",	"All2",	"Bac1",	"Cont1",	"Cont2",	"Cont3",	"Heat1",	"Heat2",	"Heat3",	"pCO2-1",	"pCO2-2",	"All1",	"All2",	"All3",	"Bac1",	"Bac2",	"Cont1",	"Cont2",	"Cont3",	"Heat1",	"Heat2",	"pCO2-1",	"pCO2-2",	"pCO2-3",	"All1",	"All2",	"All3",	"Bac1",	"Bac2",	"Bac3",	"Cont1",	"Cont2",	"Cont3",	"Heat1",	"Heat2",	"Heat3",	"pCO2-1",	"pCO2-2",	"pCO2-3")

sample_data(NewSampleNames.ps)$Genotype <- factor(sample_data(NewSampleNames.ps)$Genotype, levels= c("G4", "G12", "G20", "G38", "G2", "G27", "G31", "G34"))

Plot1AOver.1 <- plot_bar(NewSampleNames.ps, x="sample_Sample", fill="Class") + 
  facet_wrap("Genotype", scales="free_x", nrow = 2) + labs(y="Relative Abundance") + 
  theme(panel.spacing.y= unit(4, "lines"), panel.spacing = unit(2, "lines"), axis.text.x=element_text(vjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ geom_bar(stat = "identity") + scale_x_discrete(limits=c("Cont1",	"Cont2",	"Cont3","Bac1",	"Bac2",	"Bac3","Heat1",	"Heat2",	"Heat3",	"pCO2-1",	"pCO2-2",	"pCO2-3",	"All1",	"All2", "All3")) +
  scale_fill_manual(values=c("green3", "forestgreen", "darkgreen", "#00CCFF", "#0066FF", "#FFFF00"), 
                    name="Symbiont genus and type", breaks=c("C1232", "C3", "C3k", "D1", "D1a"),
                    labels=c("Cladocopium 1232", "Cladocopium 3", "Cladocopium 3k", "Durusdinium D1","Durusdinium D1a"))
Plot1AOver.1

Plot1AOver.1 <- plot_bar(NewSampleNames.ps, x="sample_Sample", fill="Class") + 
  facet_wrap("Genotype", scales="free_x", nrow = 2) + labs(y="Relative Abundance") + 
  theme(panel.spacing.y= unit(4, "lines"), panel.spacing = unit(2, "lines"), axis.text.x=element_text(vjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ geom_bar(stat = "identity") + scale_x_discrete(limits=c("Cont1",	"Cont2",	"Cont3","Bac1",	"Bac2",	"Bac3","Heat1",	"Heat2",	"Heat3",	"pCO2-1",	"pCO2-2",	"pCO2-3",	"All1",	"All2", "All3")) +
  scale_fill_manual(values=c("#66CC33", "forestgreen", "darkgreen", "#00CCFF", "#0066FF", "#FFFF00"), 
                    name="Symbiont genus and type", breaks=c("C1232", "C3", "C3k", "D1", "D1a"),
                    labels=c("Cladocopium 1232", "Cladocopium 3", "Cladocopium 3k", "Durusdinium D1","Durusdinium D1a"))
Plot1AOver.1

Plot1AOver.1 <- plot_bar(NewSampleNames.ps, x="sample_Sample", fill="Class") + 
  facet_wrap("Genotype", scales="free_x", nrow = 2) + labs(y="Relative Abundance") + 
  theme(panel.spacing.y= unit(4, "lines"), panel.spacing = unit(2, "lines"), axis.text.x=element_text(vjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ geom_bar(stat = "identity") + scale_x_discrete(limits=c("Cont1",	"Cont2",	"Cont3","Bac1",	"Bac2",	"Bac3","Heat1",	"Heat2",	"Heat3",	"pCO2-1",	"pCO2-2",	"pCO2-3",	"All1",	"All2", "All3")) +
  scale_fill_manual(values=c("#33CC33", "forestgreen", "darkgreen", "#00CCFF", "#0066FF", "#FFFF00"), 
                    name="Symbiont genus and type", breaks=c("C1232", "C3", "C3k", "D1", "D1a"),
                    labels=c("Cladocopium 1232", "Cladocopium 3", "Cladocopium 3k", "Durusdinium D1","Durusdinium D1a"))
Plot1AOver.1

###########################################################
##Normalized beta diversity metrics to compare treatments##
###########################################################
#start with ps from beginning
qd.full <- ps

#transform counts data using DESeq2
diagdds = phyloseq_to_deseq2(qd.full, ~ Genotype)

#rlog transform
rlogCOUNTS<-rlog(diagdds,blind=TRUE)
head(assay(rlogCOUNTS))

dat=as.data.frame(assay(rlogCOUNTS)) 

dat[dat < 0.0] <- 0.0 #remove negative values before calculating bray curtis distances
otu_table(qd.full) <- otu_table(dat, taxa_are_rows = TRUE) #put this transformed data back into phyloseq object as OTU table
write.csv(otu_table(qd.full),file="GCDec19_ASVs_DEseq_LONG.csv",quote=F)

# Calculate all four distance matrices: weighted unifrac, unweighted unifrac, bray curtis, binary jaccard
GC_bc <- phyloseq::distance(ps, method = "bray")
GC_bj <- phyloseq::distance(ps, method = "jaccard", binary =TRUE)

## Add between group distances into a mapping file for plotting
bc_means <- data.frame(rowMeans(as.matrix(GC_bc)))
bj_means <- data.frame(rowMeans(as.matrix(GC_bj)))
# Make a dataframe
# Extracts the SampleIDs from the dataframe
SampleID <- row.names(bc_means) 
# makes a new dataframe with the specified columnes
raw_distances <- data.frame(SampleID, bc_means,bj_means) 
# Makes metadata into a df to work with
s <- data.frame(sample_data(GC)) 
# Change first column title to SampleID to match distances dataframe
colnames(s)[1] <- "SampleID" 
# merges metadata df and distances df
bet_distances <- merge(raw_distances, s, by = "SampleID") 
colnames(bet_distances)[2:3] <- c("betdistbc","betdistbj")

write.csv(bet_distances, file = "/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/R2/bet_distances_deseq2_TWO.csv")

#I divided each individual sample bray curtis distance in above file by the average bray curtis distance of the controls within the same genotype
#now reading back in the "normalized" change in bc distance

#read in modified beta diveristy distance file
bet_mod <- read.csv(file = "/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/R2/bet_distances_deseq2_mod.csv")

##determining averages for best and worst controls
controls <- bet_mod %>%
  filter(Treatment2== "control")
tapply(controls$betdistbc, controls$Status, mean)

stress <- bet_mod %>%
  filter(Treatment2== "stress")
tapply(stress$betdistbc, stress$Status, mean)
shapiro.test(stress$betdistbc_mod) #not normal
wilcox_stat <- wilcox.test(betdistbc_mod ~ Status, data=stress)
wilcox_stat #sig difference bw best and worst normalized change in beta diversity in stress fragments

worst <- bet_mod %>%
  filter(Status== "Worst")
shapiro.test(worst$betdistbc_mod) 
pairwise.wilcox.test(worst$betdistbc_mod, worst$Treatment, p.adjust.method = "BH")

best <- bet_mod %>%
  filter(Status== "Best")
shapiro.test(best$betdistbc_mod) 
pairwise.wilcox.test(best$betdistbc_mod, best$Treatment)

bet_mod %>%
  filter(Status== "Worst") %>%
  summarise_each(funs(wilcox.test(.[Treatment == "Heat"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)

bet_mod %>%
  filter(Status== "Worst") %>%
  summarise_each(funs(wilcox.test(.[Treatment == "Bacteria"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)

##Stats for between-group distances by overall worst and best (check for both)
bet_mod %>%
  filter(Status== "Worst") %>%
  summarise_each(funs(t.test(.[Treatment == "Heat"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)

bet_mod %>%
  filter(Status== "Worst") %>%
  summarise_each(funs(t.test(.[Treatment == "All"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)

bet_mod %>%
  filter(Status== "Worst") %>%
  summarise_each(funs(t.test(.[Treatment == "Bacteria"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)

bet_mod %>%
  filter(Status== "Worst") %>%
  summarise_each(funs(t.test(.[Treatment == "pCO2"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)

bet_mod$Genotype <- factor(bet_mod$Genotype, levels= c("G4", "G12", "G20", "G38", "G2", "G27", "G31", "G34"))
bet_mod$Treatment <- factor(bet_mod$Treatment,levels = c("Control", "Bacteria", "Heat", "pCO2", "All"))

labels <- c(G4 = "Best-4", G12 ="Best-12", G20="Best-20", G38="Best-38", G2="Worst-2", G27="Worst-27", G31="Worst-31", G34="Worst-34")
p1bc <- ggplot(data = bet_mod, aes(x=Treatment, y=betdistbc_mod)) +
  geom_boxplot(alpha=0.5) +
  facet_wrap("Genotype", nrow=2, labeller= labeller(Genotype=labels)) + 
  xlab("") +
  scale_x_discrete(labels=c("Control", "Bacteria", "Heat", "pCO2", "Combined")) +
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.background = element_rect(fill="#FFFFFF")) +
  theme(panel.background = element_rect(fill="white"),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.title=element_blank(),
        axis.text = element_text(size =8),
        axis.text.x=element_text(angle = -90, vjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position="none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )
p1bc

labels <- c(Worst = "Worst Genets Average", Best = "Best Genets Average")
p2bc <- ggplot(data = bet_mod, aes(x=Treatment, y=betdistbc_mod)) +
  geom_boxplot(alpha=0.5) + 
  facet_wrap("Status", nrow=2, labeller= labeller(Status=labels)) + 
  scale_x_discrete(labels=c("Control", "Bacteria", "Heat", "pCO2", "Combined")) +
  theme(strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(fill="#FFFFFF")) +
  xlab("") + ylab("Normalized Bray Curtis between-group distance") +
  theme(panel.background = element_rect(fill="white"),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size=1),
        legend.title=element_blank(),
        axis.text = element_text(size =8),
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position="none"
  )
p2bc
library(ggpubr)
ggarrange(p2bc, p1bc, widths = c(1.25, 2), 
          labels = c("A", "B"),
          ncol = 2)

ggsave("/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/R2/allnormalizedbetas.v2.jpeg", width=10, height=6)
ggsave("/Users/laurenhowe-kerr/Dropbox/Genetic_Constraints_LHK/R2/allnormalizedbetas.v2.eps", width=10, height=6)

Worst27 <- bet_mod %>%
  filter(Rank== 8)

Worst31 <- bet_mod %>%
  filter(Rank== 7)

Worst34 <- bet_mod %>%
  filter(Rank== 6)

Worst2 <- bet_mod %>%
  filter(Rank== 5)

Best20 <- bet_mod %>%
  filter(Rank== 4)

Best12 <- bet_mod %>%
  filter(Rank== 3)

Best4 <- bet_mod %>%
  filter(Rank== 2)

Best38 <- bet_mod %>%
  filter(Rank== 1)

shapiro.test(Worst31$betdistbc_mod)
shapiro.test(Worst27$betdistbc_mod) 


# Compute the analysis of variance
res.aov <- aov(betdistbc_mod ~ Treatment, data = Worst31)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

pairwise.t.test(Worst31$betdistbc_mod, Worst31$Treatment)
pairwise.t.test(Best38$betdistbc_mod, Best38$Treatment)
 
#test differences by genet
bet_mod %>%
  filter(Rank== 1) %>%
  summarise_each(funs(t.test(.[Treatment == "All"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)
bet_mod %>%
  filter(Rank== 1) %>%
  filter(Treatment=="All" & Treatment== "Control")

bet_mod %>%
  filter(Rank== 1) %>%
  summarise_each(funs(t.test(.[Treatment == "Control"], .[Treatment == "Heat"])$p.value), vars = betdistbc_mod)
bet_mod %>%
  filter(Rank== 1) %>%
  summarise_each(funs(t.test(.[Treatment == "Bacteria"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)
bet_mod %>%
  filter(Rank== 1) %>%
  summarise_each(funs(t.test(.[Treatment == "pCO2"], .[Treatment == "Control"])$p.value), vars = betdistbc_mod)

                             
###########################################################
##########################GJAM analyses####################
###########################################################      
                             
 x<-read.csv("GJAM_GC19_sampledata_LULU.csv",header=TRUE) # Load the covariates
y<-read.csv("GJAM_GC19_OTUs_LULU.csv",header=TRUE) #Load the OTU table

# Make sure y and x samples are organized in the same order
y[,1]<-as.character(y[,1])
x[,1]<-as.character(x[,1])
y<-y[order(y[,1]),]
x<-x[order(x[,1]),]

#rearrange the treatment and genotype factors to have the right reference (for example control treatment as the reference treatment)
x[,4]<-as.character(x[,4])  
x[,4][x[,4]=="All"]<-"Combined"
x[,4]<-as.factor(x[,4])
x[,4]<-relevel(x[,4],"Control")
x[,4]<-factor(x[,4],levels=c("Control","Combined","Bacteria","Heat","pCO2"))
x[,3]<-as.character(x[,3])
x[,3][x[,3]=="G38"]<-"Best-38"
x[,3][x[,3]=="G12"]<-"Best-12"
x[,3][x[,3]=="G2"]<-"Worst-12"
x[,3][x[,3]=="G20"]<-"Best-20"
x[,3][x[,3]=="Worst-12"]<-"Worst-2"
x[,3][x[,3]=="G27"]<-"Worst-27"
x[,3][x[,3]=="G31"]<-"Worst-31"
x[,3][x[,3]=="G34"]<-"Worst-34"
x[,3][x[,3]=="G4"]<-"Best-4"
x[,3]<-as.factor(x[,3])
x[,3]<-factor(x[,3],levels=c("Best-38","Best-4","Best-12","Best-20","Worst-2","Worst-34”,”Worst-31”,”Worst-27"))
x[,3]<-relevel(x[,3],ref="Best-38")

# remove the sample column
y<-y[,-1]

# Load GJAM
library(gjam)

#Parameter for the models (number of iterations, length of burnin, and types of data)
ml   <- list(ng = 10000, burnin = 500, typeNames = "CC")

# Run the regional model (across reefs)
out  <- gjam(~Treatment+Genotype, x, y, modelList = ml)
full <- gjamSensitivity(out)

# Run local models (within reef)
yPandora <-y[x$Genotype%in%c("Worst-27","Worst-31","Worst-34"),]
xPandora <-x[x$Genotype%in%c("Worst-27","Worst-31","Worst-34"),]
xPandora[,3]<-as.factor(as.character(xPandora[,3]))
xPandora[,3]<-factor(xPandora[,3],levels=c("Worst-34","Worst-31","Worst-27"))
outPandora  <- gjam(~Treatment+Genotype, xPandora, yPandora, modelList = ml)
Pandora <- gjamSensitivity(outPandora)

yRib <-y[x$Genotype%in%c("Best-38","Best-12","Best-20"),]
xRib <-x[x$Genotype%in%c("Best-38","Best-12","Best-20"),]
xRib[,3]<-as.factor(as.character(xRib[,3]))
xRib[,3]<-factor(xRib[,3],levels=c("Best-38","Best-12","Best-20"))
outRib  <- gjam(~Treatment+Genotype, xRib, yRib, modelList = ml)
Rib <- gjamSensitivity(outRib)

