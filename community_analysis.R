library(tidyverse)
library(vegan)
library(reshape2)
library(ggrepel)


# Set this to wherever github repo is saved
path<-"./scripts_repo"
setwd(path)
list.files(".")

# Source functions taken from here:
#https://oliviarata.wordpress.com/2014/07/17/ordinations-in-ggplot2-v2-ordisurf/

source("./ord_labels.R")
source("./gg_ordisurf.R")

# Read in taxa abundances:
df_mat<-read.csv("./metaxa2_output_matrix.csv", row.names = 1)
head(df_mat)

# Check relative abundances sums to 1:
min(rowSums(df_mat))
max(rowSums(df_mat))

# Get metadata:

df2<-read.csv("./data.csv")
head(df2)

# Reorder metadata to match matrix. Merge, split, check.
tmp<-merge(df2, df_mat, by.x = "sample", by.y = "row.names")
dim(tmp)
row.names(tmp)<-tmp$sample
head(tmp)
mat<-tmp[,10:116]
meta<-tmp[,1:9]

all.equal(row.names(mat), row.names(meta))

# Permutational anova
# Variance in DR abundance explained by class-level taxa abundances.
# Patchy distribution at lower phylo levels.
class_res<-adonis(mat ~ DR_count, data = meta, permutations = 9999, method = "bray") 
class_res


# Plot community overview. How do they cluster? Which taxa drive each axis?

all.mds <- metaMDS(mat, k = 4, distance = "bray")

# Get species loadings:
# Extract top loadings...
NMDS = data.frame(MDS1 = all.mds$points[,1], MDS2 = all.mds$points[,2])
vec.sp<-envfit(all.mds$points, mat, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)
head(vec.sp.df)

# Get just top 10 species:
tmp<-merge(as.data.frame(vec.sp$vectors$r), as.data.frame(vec.sp$vectors$pvals), by = "row.names")
colnames(tmp)<-c("species", "r2", "pval")
head(tmp)

# Drop uninformative species:
drops<-c('unassigned', 'Family Incertae Sedis')
top10<-tmp %>%
  subset(., !species %in% drops) %>%
  arrange(., -r2) %>%
  head(n = 10) %>%
  select(., species)
top10<-unique(top10$species)

vec.sp.df.10<-subset(vec.sp.df, species %in% top10)

NMDS %>%
  mutate(., sample = row.names(.)) %>%
  left_join(., meta, by = "sample") %>%
  ggplot(., aes(MDS1, MDS2)) + 
  geom_point(size = 1, alpha = 0.3, color = "grey2")+
  geom_segment(data=vec.sp.df.10,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")), color = "grey")+ 
  geom_text_repel(data=vec.sp.df.10,aes(x=MDS1,y=MDS2,label=species),size=3, color = "black")+
  coord_fixed(ratio=1)+
  theme_classic()+
  xlab("NMDS1")+
  ylab("NMDS2")+
  theme(text = element_text(size = 15),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

# Or with EMPO levels overlaid:
comm_plot<-NMDS %>%
  mutate(., sample = row.names(.)) %>%
  left_join(., meta, by = "sample") %>%
  ggplot(., aes(MDS1, MDS2)) + 
  geom_point(size = 1, alpha = 0.8, aes(color = EMPO_level2))+
  geom_segment(data=vec.sp.df.10,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")), color = "grey")+ 
  geom_text_repel(data=vec.sp.df.10,aes(x=MDS1,y=MDS2,label=species),size=3, color = "black")+
  coord_fixed(ratio=1)+
  theme_classic()+
  xlab("NMDS1")+
  ylab("NMDS2")+
  theme(text = element_text(size = 15))+
        #axis.title.x=element_blank(),
        #axis.title.y=element_blank())+
  scale_color_brewer(type = "qual", palette = 7)


# So are certain communities enriched for CRISPR?
# Plot the results:

head(all.mds)
stressplot(all.mds)

a<-plot(all.mds)
head(a)
sample_data<-as.data.frame(a$sites)
sample_data$sample<-row.names(sample_data)

# Check CRISPR abundance (GAM model preidctions) and overlay on MDS 

dr_plot<-gg_ordisurf(all.mds, env.var = meta$DR_count, pt.size = 1, binwidth = 15, var.label = "CRISPR\nAbundance")

head(dr_plot$plot)
dr_plot_output<-dr_plot$plot
# Panel plot:

egg::ggarrange(comm_plot, dr_plot_output, ncol = 2, labels = c('A', 'B'))

