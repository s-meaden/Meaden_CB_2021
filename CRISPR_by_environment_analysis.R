library(tidyverse)
library(reshape2)
library(broom)

# Descriptive statistics of CRISPR abundance variation across environments.

# Set this to wherever github repo is saved
path<-"./scripts_repo"

full_df<-read.csv("./data.csv")
head(full_df)

# Variation explained by EMPO level 1
hist(log10(full_df$DR_count))
m1<-glm(log10(DR_count) ~ EMPO_level1, data = full_df)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
m2<-update(m1,~.- EMPO_level1)
anova(m1, m2, test = "F")
lev1<-anova(m1, m2, test = "F")
with(summary(m1), 1 - deviance/null.deviance)

# Host associated environments have significantly more CRISPR arrays than free-living.


m1<-glm(log10(DR_count) ~ EMPO_level2, data = full_df)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
m2<-update(m1,~.- EMPO_level2)
anova(m1, m2, test = "F")
lev2<-anova(m1, m2, test = "F")
# Adjust p-val:
p.adjust(0.001291, method = 'bonferroni', n = 3) 
with(summary(m1), 1 - deviance/null.deviance)

# Level 3 variation

m1<-glm(log10(DR_count) ~ EMPO_level3, data = full_df)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
m2<-update(m1,~.- EMPO_level3)
anova(m1, m2, test = "F")
lev3<-anova(m1, m2, test = "F")
p.adjust(1.596e-14, method = 'bonferroni', n = 3)
with(summary(m1), 1 - deviance/null.deviance)


# Table of results:
lev1<-lev1 %>%
  tidy() %>%
  mutate(., EMPO_level = 1) %>%
  slice(2)
lev2<-lev2 %>%
  tidy() %>%
  mutate(., EMPO_level = 2) %>%
  slice(2)
lev3<-lev3 %>%
  tidy() %>%
  mutate(., EMPO_level = 3) %>%
  slice(2)

tmp<-lev1 %>%
  bind_rows( lev2, lev3) %>%
  select( Resid..Df, statistic, p.value, EMPO_level) %>%
  rename( F_value = statistic, DF = Resid..Df) %>%
  mutate( corrected_p_value = p.adjust(p.value, method = 'bonferroni', n = 3)) %>%
  mutate( F_value = round( F_value, digits = 4),
          p.value = round( p.value, digits = 4),
          corrected_p_value = round(corrected_p_value, digits = 4)) %>%
  mutate( p.value = case_when(p.value == "0" ~ "<0.0001",
                              TRUE ~ as.character(p.value)),
          corrected_p_value = case_when(corrected_p_value == "0" ~ "<0.0001",
                                        TRUE ~ as.character(corrected_p_value)))

# Panel plot of results:
a<-full_df %>%
  select(sample, DR_count, EMPO_level1) %>%
  melt(., id.vars  = c('sample', 'DR_count')) %>%
  ggplot(., aes(value, DR_count))+
  geom_jitter(width = 0.05, size = 1)+
  geom_boxplot( alpha = 0.6)+
  scale_y_log10()+
  xlab("")+
  ylab("CRISPR\nAbundance")+
  theme_classic()+
  ggtitle("EMPO Level 1")+
  theme(text = element_text(size = 15))

b<-full_df %>%
  select(sample, DR_count, EMPO_level2) %>%
  melt(., id.vars  = c('sample', 'DR_count')) %>%
  ggplot(., aes(value, DR_count))+
  geom_jitter(width = 0.05, size = 1)+
  geom_boxplot( alpha = 0.6)+
  scale_y_log10()+
  xlab("")+
  ylab("CRISPR\nAbundance")+
  theme_classic()+
  ggtitle("EMPO Level 2")+
  theme(text = element_text(size = 15))

# Drop envs with < 3 samples.
full_df %>%
  select(sample, DR_count, EMPO_level3) %>%
  melt(., id.vars  = c('sample', 'DR_count')) %>%
  group_by(value) %>%
  tally() %>%
  arrange(n) %>%
  head()
drops<-c('animal_secretion', 'proximal_gut')

# Order by descending abundance
tmp<-full_df %>%
  select(sample, DR_count, EMPO_level3) %>%
  filter(., !EMPO_level3 %in% drops) %>%
  melt(., id.vars  = c('sample', 'DR_count')) %>%
  group_by(value) %>%
  summarise(mean = mean(DR_count)) %>%
  arrange( -mean)
positions<-tmp$value

positions<-gsub("_", "\n", positions)

c<-full_df %>%
  select(sample, DR_count, EMPO_level3) %>%
  filter(., !EMPO_level3 %in% drops) %>%
  melt(., id.vars  = c('sample', 'DR_count')) %>%
  mutate(., value = gsub("_", "\n", value)) %>%
  ggplot(., aes(value, DR_count))+
  geom_jitter(width = 0.05, size = 1)+
  geom_boxplot( alpha = 0.6)+
  scale_y_log10()+
  xlab("Environment")+
  ylab("CRISPR\nAbundance")+
  theme_classic()+
  ggtitle("EMPO Level 3")+
  theme(text = element_text(size = 15))+
  scale_x_discrete(limits = positions)
c
cowplot::plot_grid(a, b, c, ncol = 1)

## Updated figure for MS. 
# Just distributions.

a<-full_df %>%
  group_by(., EMPO_level3) %>%
  mutate(., mean = mean(DR_count)) %>%
  ungroup() %>%
  select(., sample, DR_count, EMPO_level3, mean) %>%
  mutate(EMPO_level3 = fct_reorder(EMPO_level3, desc(mean))) %>%
  filter(., !EMPO_level3 %in% drops) %>%
  ggplot(., aes(DR_count, EMPO_level3))+
  ggridges::geom_density_ridges(scale = 2, aes(fill = EMPO_level3), alpha = 0.6,
                                jittered_points = TRUE, point_size = 0.3)+
  theme_classic()+
  xlab("CRISPR Abundance\n(read count)")+
  ylab("Environment type\n(EMPO Level 3)")+
  guides(fill = FALSE)+
  #scale_fill_viridis_d(option = "magma")
  #scale_fill_viridis_d()+
  scale_fill_brewer(type = "qual", palette = 3)+
  theme(text = element_text(size = 15))

b<-full_df %>%
  group_by(., EMPO_level2) %>%
  mutate(., mean = mean(DR_count)) %>%
  ungroup() %>%
  select(., sample, DR_count, EMPO_level2, mean) %>%
  mutate(EMPO_level2 = fct_reorder(EMPO_level2, desc(mean))) %>%
  ggplot(., aes(DR_count, EMPO_level2))+
  ggridges::geom_density_ridges(scale = 2, aes(fill = EMPO_level2), alpha = 0.6,
                                jittered_points = TRUE, point_size = 0.3)+
  theme_classic()+
  xlab("CRISPR Abundance\n(read count)")+
  ylab("Environment type\n(EMPO Level 2)")+
  guides(fill = FALSE)+
  #scale_fill_viridis_d(option = "magma")
  #scale_fill_brewer(type = "div", palette = 1)+
  scale_fill_brewer(type = "qual", palette = 1)+
  theme(text = element_text(size = 15))

c<-full_df %>%
  group_by(., EMPO_level1) %>%
  mutate(., mean = mean(DR_count)) %>%
  ungroup() %>%
  select(., sample, DR_count, EMPO_level1, mean) %>%
  mutate(EMPO_level1 = fct_reorder(EMPO_level1, desc(mean))) %>%
  ggplot(., aes(DR_count, EMPO_level1))+
  ggridges::geom_density_ridges(scale = 5, aes(fill = EMPO_level1), alpha = 0.6,
                                jittered_points = TRUE, point_size = 0.3)+
  theme_classic()+
  xlab("CRISPR Abundance\n(read count)")+
  ylab("Environment type\n(EMPO Level 1)")+
  guides(fill = FALSE)+
  #scale_fill_viridis_d(option = "magma")
  scale_fill_brewer(type = "qual", palette = 2)+
  theme(text = element_text(size = 15))

# And add placeholder panel where EMPO flowchart will go
egg::ggarrange(c, c, b, a, labels = c('A', 'B', 'C', 'D'), ncol = 1)


# Post-hoc testing of level 3:
tmp<-full_df %>%
  select(sample, DR_count, EMPO_level3) %>%
  filter(., !EMPO_level3 %in% drops)

m1<-glm(log10(DR_count) ~ EMPO_level3, data = tmp)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
m2<-update(m1,~.-EMPO_level3)
anova(m1, m2, test = "F")

##### Post-hoc testting for SI.

emm <- emmeans::emmeans(m1, "EMPO_level3", transform = "response")
emm
tmp<-tidy(pairs(emm))   # "P value adjustment: tukey method for comparing a family of 10 estimates"

tmp<-tmp %>%
  mutate(estimate = round( estimate, digits = 2),
         std.error = round( std.error, digits = 2),
         df = round( df , digits = 2),
         z.ratio = round( z.ratio, digits = 2),
         adj.p.value = round( adj.p.value, digits = 2)) %>%
  mutate( adj.p.value = case_when(adj.p.value == "0" ~ "<0.0001",
                              TRUE ~ as.character(adj.p.value)))









