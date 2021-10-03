library(tidyverse)
library(plotly)
library(sjPlot)

# Does CRISPR abundance correlate with phage diversity. 
# Theory predicts CRISPR only useful when viral div relatively low.

# Set this to wherever github repo is saved
path<-"./scripts_repo"

# Get data:
vir_load<-read.csv("./vir_abund_data.csv")

head(vir_load)

full_df<-read.csv("./data.csv")
head(full_df)

df<-full_df %>%
  left_join(., vir_load, by = "sample") %>%
  filter(., complete.cases(.))

# Add viral contig counts:
vir_div<-read.csv("./vir_nctg.csv")

# Link up
df<-df %>%
  left_join(vir_div, by = "sample")

# Check if viral load ~ viral div i.e. sampling effort.

df %>%
  ggplot(., aes(viral_load, n_ctg))+
  geom_point( color = "darkgrey", size = 1)+
  scale_x_log10(breaks = c(1e+05, 1e+02, 1e+03, 1e+04))+
  scale_y_log10()+
  geom_smooth( method = "lm", color = "black", linetype = "dashed", alpha = 0.5)+
  theme_classic()+
  xlab("Viral load\n(reads per million)")+
  ylab("Viral contigs\nper sample")

# Get R2:
hist(log10(df$viral_load))
hist(log10(df$n_ctg))

m1<-glm(log10(viral_load) ~ log10(n_ctg), data = df)
summary(m1)
with(summary(m1), 1 - deviance/null.deviance)


# Normalise diversity metrics by viral load to account for sampling effort.

# Make plot first. Then do stats:

tmp<-df %>%
  mutate(., nei_norm = nei_div / viral_load) %>%
  mutate(., rich_norm = n_ctg / viral_load) %>%
  mutate(., even_norm = evenness / viral_load) %>%
  mutate(., shan_norm = shannon_div / viral_load) %>%
  select( DR_count, nei_norm, rich_norm, even_norm, shan_norm)

tmp %>%
  rename("Nei\nIndex" = nei_norm, "Richness" = rich_norm, "Evenness" = even_norm, "Shannon\nIndex" = shan_norm) %>%
  reshape2::melt(., id.vars = "DR_count") %>%
  ggplot(., aes(value, DR_count))+
  geom_point( color = "darkgray", size = 1)+
  geom_smooth(method = "lm", linetype = "dashed", color = "black")+
  facet_wrap( ~ variable, scales = "free")+
  scale_y_log10()+
  scale_x_log10()+
  ylab("CRISPR Abundance\n(DR Count)")+
  xlab("Normalised Diversity Score")+
  theme_classic()+
  theme(text = element_text(size = 15))

# Plot for MS with panels labelled a , b , c etc:

tmp2<-tmp %>%
  rename("Nei\nIndex" = nei_norm, "Richness" = rich_norm, "Evenness" = even_norm, "Shannon\nIndex" = shan_norm) %>%
  reshape2::melt(., id.vars = "DR_count")
  
a<-tmp2 %>%
  filter( variable == "Nei\nIndex") %>%
  ggplot(., aes(value, DR_count))+
  geom_point( color = "darkgray", size = 1)+
  geom_smooth(method = "lm", linetype = "dashed", color = "black")+
  facet_wrap( ~ variable, scales = "free")+
  scale_y_log10()+
  scale_x_log10()+
  ylab("CRISPR Abundance\n(DR Count)")+
  xlab("Normalised Diversity Score")+
  theme_classic()+
  theme(text = element_text(size = 15))

b<-tmp2 %>%
  filter( variable == "Shannon\nIndex") %>%
  ggplot(., aes(value, DR_count))+
  geom_point( color = "darkgray", size = 1)+
  geom_smooth(method = "lm", linetype = "dashed", color = "black")+
  facet_wrap( ~ variable, scales = "free")+
  scale_y_log10()+
  scale_x_log10()+
  ylab("CRISPR Abundance\n(DR Count)")+
  xlab("Normalised Diversity Score")+
  theme_classic()+
  theme(text = element_text(size = 15))

c<-tmp2 %>%
  filter( variable == "Richness") %>%
  ggplot(., aes(value, DR_count))+
  geom_point( color = "darkgray", size = 1)+
  geom_smooth(method = "lm", linetype = "dashed", color = "black")+
  facet_wrap( ~ variable, scales = "free")+
  scale_y_log10()+
  scale_x_log10()+
  ylab("CRISPR Abundance\n(DR Count)")+
  xlab("Normalised Diversity Score")+
  theme_classic()+
  theme(text = element_text(size = 15))

d<-tmp2 %>%
  filter( variable == "Evenness") %>%
  ggplot(., aes(value, DR_count))+
  geom_point( color = "darkgray", size = 1)+
  geom_smooth(method = "lm", linetype = "dashed", color = "black")+
  facet_wrap( ~ variable, scales = "free")+
  scale_y_log10()+
  scale_x_log10()+
  ylab("CRISPR Abundance\n(DR Count)")+
  xlab("Normalised Diversity Score")+
  theme_classic()+
  theme(text = element_text(size = 15))

egg::ggarrange(a, b, c, d, labels = c('A', 'B', 'C', 'D'))

# Make df with results:
stats_res<-data.frame(metric = c('shannon', 'richness', 'evenness', 'nei'))

hist(log10(tmp$DR_count))
m1<-lm(log10(DR_count) ~ log10(shan_norm), data = tmp)
summary(m1)
m2<-update(m1,~.-log10(shan_norm))
anova(m1, m2, test = "F")
f_val<-anova(m1, m2, test = "F")[2,5]
p_val<-anova(m1, m2, test = "F")[2,6]

hist(log10(tmp$DR_count))
m1<-lm(log10(DR_count) ~ log10(rich_norm), data = tmp)
summary(m1)
m2<-update(m1,~.-log10(rich_norm))
anova(m1, m2, test = "F")
f_val <- c(f_val, anova(m1, m2, test = "F")[2,5])
p_val <- c(p_val, anova(m1, m2, test = "F")[2,6])

m1<-lm(log10(DR_count) ~ log10(even_norm), data = tmp)
summary(m1)
m2<-update(m1,~.-log10(even_norm))
anova(m1, m2, test = "F")
f_val <- c(f_val, anova(m1, m2, test = "F")[2,5])
p_val <- c(p_val, anova(m1, m2, test = "F")[2,6])

m1<-lm(log10(DR_count) ~ log10(nei_norm), data = tmp)
summary(m1)
m2<-update(m1,~.-log10(nei_norm))
anova(m1, m2, test = "F")
f_val <- c(f_val, anova(m1, m2, test = "F")[2,5])
p_val <- c(p_val, anova(m1, m2, test = "F")[2,6])

stats_res$F_vals<-f_val
stats_res$P_vals<-p_val

stats.res<-stats_res %>%
  mutate( p.adj = p.adjust(P_vals, method = "bonferroni", n = 4)) %>%
  mutate( P_vals = round( P_vals , digits = 3),
          F_vals = round(F_vals, digits = 3),
          p.adj = round(p.adj, digits = 3)) %>%
  mutate(p.adj = case_when(p.adj == "0" ~ "<0.001",
                          TRUE ~ as.character(p.adj)),
         P_vals = case_when(P_vals == "0" ~ "<0.001",
                   TRUE ~ as.character(P_vals))) %>%
  rename(F_value = F_vals, P_value = P_vals, Adjusted_p_value = p.adj)






