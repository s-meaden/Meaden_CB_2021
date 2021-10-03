
library(tidyverse)
library(plotly)

# Get viral load data:

# Set this to wherever github repo is saved
path<-"./scripts_repo"

vir_load<-read.csv("./vir_abund_data.csv")
head(vir_load)

# Link up with other data
full_df<-read.csv("./data.csv")

# Add in contig counts:
contigs<-read.csv("./vir_nctg.csv")

# Link all together

df<-full_df %>%
  left_join(., contigs, by = "sample") %>%
  left_join(., vir_load, by = "sample") %>%
  filter(., complete.cases(.)) %>%
  rename( viral_abund = viral_load)

# Does phage abundance predict CRISPR abundance across all environments

all<-df %>%
  ggplot(., aes(DR_count, viral_abund))+
  geom_point()+
  scale_y_log10()+
  scale_x_log10()+
  geom_smooth( method = "lm", linetype = "dashed")+
  theme_classic()+
  xlab("CRISPR abundance\n(DR count)")+
  ylab("Viral Abundance\n(sum of viral contig coverage)")+
  theme(text = element_text(size = 15))

m1<-glm(log10(DR_count) ~ log10(viral_abund), data = df)
qqnorm(resid(m1))
qqline(resid(m1))
m2<-update(m1,~.-log10(viral_abund))
anova(m1, m2, test = "F")
summary(m1)
with(summary(m1), 1 - deviance/null.deviance) 


# Plot the correlations split  by environment for supplemental / panel plot:
# Add colors manually for consistency later

pal<-c('#1b9e77', '#d95f02')
a<-df %>%
  mutate(EMPO_level1 = gsub("_", "\n", EMPO_level1)) %>%
  ggplot(., aes(DR_count, viral_abund))+
  geom_point( aes(color = EMPO_level1))+
  scale_y_log10()+
  scale_x_log10()+
  geom_smooth( method = "lm", linetype = "dashed", aes(color = EMPO_level1))+
  theme_classic()+
  xlab("CRISPR abundance\n(DR count)")+
  ylab("Viral Relative Abundance")+
  theme(text = element_text(size = 15))+
  labs(color = "EMPO\nLevel 1")+
  scale_color_manual(values = pal)
  scale_color_brewer(type = "qual", palette = 2)


pal<-c('#e41a1c', '#377eb8', '#4daf4a')
b<-df %>%
  mutate(EMPO_level2 = gsub("_", "\n", EMPO_level2)) %>%
  ggplot(., aes(DR_count, viral_abund))+
  geom_point( aes(color = EMPO_level2))+
  scale_y_log10()+
  scale_x_log10()+
  geom_smooth( method = "lm", linetype = "dashed", aes(color = EMPO_level2))+
  theme_classic()+
  xlab("CRISPR abundance\n(DR count)")+
  ylab("Viral Relative Abundance")+
  theme(text = element_text(size = 15))+
  labs(color = "EMPO\nLevel 2")+
  scale_color_manual(values = pal)
  scale_color_brewer(type = "qual", palette = 6)

### Fix palette to match with the other plots

cats<-unique(df$EMPO_level3)
  
cats<-cats[!is.na(cats)]
  
col_df<-data.frame(color = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                               '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                               '#FFDB6D'), EMPO_level3 = cats)
  
col <- as.character(col_df$color)
names(col) <- as.character(col_df$EMPO_level3)
names(col)<-gsub("_", "\n", names(col))

# Split level 3 into significant positive and non-cors.

posi<-c('distal\ngut', 'water\nsaline', 'surface\nsaline')

options(scipen=0)

c<-df %>%
  mutate(EMPO_level3 = gsub("_", "\n", EMPO_level3)) %>%
  mutate( significant = ifelse(EMPO_level3 %in% posi, "adjusted p < 0.05", "NS")) %>%
  ggplot(., aes(DR_count, viral_abund))+
  geom_point( aes(color = EMPO_level3))+
  scale_y_log10()+
  scale_x_log10()+
  geom_smooth( method = "lm", linetype = "dashed", aes(color = EMPO_level3), alpha = 0.3)+
  theme_classic()+
  xlab("CRISPR abundance\n(DR count)")+
  ylab("Viral Relative Abundance")+
  theme(text = element_text(size = 15))+
  labs(color = "EMPO\nLevel 3")+
  scale_color_manual(values = col)+
  facet_wrap( ~ significant)

# Tiy up in Affinity Designer

egg::ggarrange(a, b, labels = c('A','B'), ncol= 2)
egg::ggarrange(c, labels = 'C')


