
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

# Run GLM on each environmental level- consistency across levels?

# Level 1
df_glm_1 <- df %>%
  group_by(EMPO_level1) %>%
  do(mod = glm(log10(DR_count) ~ log10(viral_abund),data = .))

# Get coefficients from list:
df_coef_1 <- df_glm_1 %>%
  do(data.frame(
    EMPO_level = .$EMPO_level1,
    var = names(coef(.$mod)),
    coef(summary(.$mod)),
    level = 1)
  )

# Repeat for EMPO level 2 and 3

df_glm_2 <- df %>%
  group_by(EMPO_level2) %>%
  do(mod = glm(log10(DR_count) ~ log10(viral_abund),data = .))

df_coef_2 <- df_glm_2 %>%
  do(data.frame(
    EMPO_level = .$EMPO_level2,
    var = names(coef(.$mod)),
    coef(summary(.$mod)),
    level = 2)
  )

df_glm_3 <- df %>%
  group_by(EMPO_level3) %>%
  do(mod = glm(log10(DR_count) ~ log10(viral_abund),data = .))

df_coef_3 <- df_glm_3 %>%
  do(data.frame(
    EMPO_level = .$EMPO_level3,
    var = names(coef(.$mod)),
    coef(summary(.$mod)),
    level = 3)
  )

# Stick them together, remove intercepts, FDR correct, label
tmp<-df_coef_1 %>%
  bind_rows(df_coef_2) %>%
  bind_rows(df_coef_3) %>%
  filter( !var == "(Intercept)")

tmp$padj<-p.adjust(tmp$Pr...t.., method = "fdr", n = length(tmp$Pr...t..))

tmp<-tmp %>%
  mutate(sig = case_when(padj <= 0.05 & padj > 0.01 ~ "*",
                         padj <= 0.01 & padj > 0.001 ~ "**",
                         padj <= 0.001 ~ "***",
                         padj > 0.05 ~ "",
                         TRUE ~ "NA")) %>%
  select( -var)

colnames(tmp)<-c("EMPO category","Estimate", "Std. Error", "t-value", "P-value", "EMPO level", "padj", "sig")

tmp<-tmp %>%
  mutate(Estimate = round( Estimate, digits = 3),
         `Std. Error` = round(`Std. Error`, 3),
         `t-value` = round( `t-value`, 3),
         `P-value`= round(`P-value`, 3),
         padj = round(padj, 3)) %>%
  mutate( padj = case_when(padj == "0" ~ "<0.001",
                           TRUE ~ as.character(padj))) %>%
  mutate( `P-value` = case_when(`P-value` == "0" ~ "<0.001",
                                TRUE ~ as.character(`P-value`))) %>%
  rename( 'adjusted-p-value' = padj ,
          significance = sig)

tmp

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

# Split level 3 into significant and non-significant correlations from post-hoc testing

posi<-c('distal\ngut', 'water\nsaline', 'surface\nsaline', 'sediment\nnonsaline')

options(scipen=0)

c<-df %>%
  mutate(EMPO_level3 = gsub("_", "\n", EMPO_level3)) %>%
  mutate( significant = ifelse(EMPO_level3 %in% posi, "adjusted p < 0.05", "NS")) %>%
  ggplot(., aes(DR_count, viral_abund, label = sample))+
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

library(plotly)
plotly::ggplotly(c)

