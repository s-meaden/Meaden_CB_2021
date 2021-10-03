library(tidyverse)
library(reshape2)

# Does viral abundance and diversity vary by environment. 
# Basic descriptive statistics.

# Set this to wherever github repo is saved
path<-"PATH_TO/scripts_repo"

vir_load<-read.csv("./vir_abund_data.csv")
head(vir_load)

full_df<-read.csv("./data.csv")
head(full_df)

# Link viral load with metadata:
df<-full_df %>%
  left_join(., vir_load, by = "sample") %>%
  rename(viral_abund = viral_load)

# Drop NA categories
df<-df %>%
  filter(., !is.na(EMPO_level3))

# Greater variation between environments than between?

# Level 1 first.
hist(log10(df$viral_abund))
m1<-glm(log10(viral_abund) ~ EMPO_level1, data = df)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
m2<-update(m1,~.- EMPO_level1)
anova(m1, m2, test = "F")
with(summary(m1), 1 - deviance/null.deviance)

# Level 2 next:
m1<-glm(log10(viral_abund) ~ EMPO_level2, data = df)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
m2<-update(m1,~.- EMPO_level2)
anova(m1, m2, test = "F")
with(summary(m1), 1 - deviance/null.deviance)

# Level 3:
m1<-glm(log10(viral_abund) ~ EMPO_level3, data = df)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
m2<-update(m1,~.- EMPO_level3)
anova(m1, m2, test = "F")
with(summary(m1), 1 - deviance/null.deviance)

# Panel plot for each EMPO level
a<-df %>%
  select(sample, viral_abund, EMPO_level1) %>%
  melt(., id.vars  = c('sample', 'viral_abund')) %>%
  ggplot(., aes(value, viral_abund))+
  geom_jitter(width = 0.05)+
  geom_boxplot( alpha = 0.6)+
  scale_y_log10()+
  xlab("")+
  ylab("Viral Load")+
  theme_classic()+
  ggtitle("EMPO Level 1")+
  theme(text = element_text(size = 15))

b<-df %>%
  select(sample, viral_abund, EMPO_level2) %>%
  melt(., id.vars  = c('sample', 'viral_abund')) %>%
  ggplot(., aes(value, viral_abund))+
  geom_jitter(width = 0.05)+
  geom_boxplot( alpha = 0.6)+
  scale_y_log10()+
  xlab("")+
  ylab("Viral Load")+
  theme_classic()+
  ggtitle("EMPO Level 2")+
  theme(text = element_text(size = 15))

# Drop envs with few samples for plots:
df %>%
  select(sample, viral_abund, EMPO_level3) %>%
  melt(., id.vars  = c('sample', 'viral_abund')) %>%
  group_by(value) %>%
  tally() %>%
  arrange(n) %>%
  head()
drops<-c('animal_secretion', 'proximal_gut')


# Plot as ggridges

df<-df %>%
  mutate(EMPO_level1 = gsub("_", "\n", EMPO_level1),
         EMPO_level2 = gsub("_", "\n", EMPO_level2),
         EMPO_level3 = gsub("_", "\n", EMPO_level3))

drops<-c('animal\nsecretion', 'proximal\ngut')

tmp<-df %>%
  select(sample, viral_abund, EMPO_level1) %>%
  filter(., !EMPO_level1 %in% drops) %>%
  melt(., id.vars  = c('sample', 'viral_abund')) %>%
  group_by(value) %>%
  summarise(mean = median(viral_abund, na.rm = TRUE)) %>%
  arrange( -mean)

tmp
positions<-tmp$value


a1<-df %>%
  select(sample, viral_abund, EMPO_level1) %>%
  ggplot(., aes(viral_abund, EMPO_level1))+
  ggridges::geom_density_ridges( scale = 5, alpha = 0.7, jittered_points = TRUE,
                                 point_size = 0.3, aes(fill = EMPO_level1))+
  scale_fill_brewer(type = "seq", palette = 7, direction = -1)+
  theme_classic()+
  xlab("Viral Load (sum of viral contig coverage depth)")+
  ylab("Density")+
  scale_y_discrete(limits = positions)+
  labs(fill = "EMPO\nLevel 1")+
  guides(fill = FALSE)

tmp<-df %>%
  select(sample, viral_abund, EMPO_level2) %>%
  filter(., !EMPO_level2 %in% drops) %>%
  melt(., id.vars  = c('sample', 'viral_abund')) %>%
  group_by(value) %>%
  summarise(mean = median(viral_abund, na.rm = TRUE)) %>%
  arrange( -mean)
positions<-tmp$value

b1<-df %>%
  select(sample, viral_abund, EMPO_level2) %>%
  ggplot(., aes(viral_abund, EMPO_level2))+
  ggridges::geom_density_ridges( scale = 5, alpha = 0.7, jittered_points = TRUE,
                                 point_size = 0.3, aes(fill = EMPO_level2))+
  scale_fill_brewer(type = "seq", palette = 2)+
  theme_classic()+
  xlab("Viral Load (sum of viral contig coverage depth)")+
  ylab("Density")+
  scale_y_discrete(limits = positions)+
  labs(fill = "EMPO\nLevel 2")+
  guides(fill = FALSE)

tmp<-df %>%
  select(sample, viral_abund, EMPO_level3) %>%
  filter(., !EMPO_level3 %in% drops) %>%
  melt(., id.vars  = c('sample', 'viral_abund')) %>%
  group_by(value) %>%
  summarise(mean = median(viral_abund, na.rm = TRUE)) %>%
  arrange( -mean)
positions<-tmp$value

c1<-df %>%
  select(sample, viral_abund, EMPO_level3) %>%
  filter(., !EMPO_level3 %in% drops) %>%
  ggplot(., aes(viral_abund, EMPO_level3))+
  ggridges::geom_density_ridges( scale = 5, alpha = 0.7, jittered_points = TRUE,
                                 point_size = 0.3, aes(fill = EMPO_level3))+
  scale_fill_brewer(type = "seq", palette = 1)+
  theme_classic()+
  xlab("Viral Load (sum of viral contig coverage depth)")+
  ylab("Density")+
  scale_y_discrete(limits = positions)+
  labs(fill = "EMPO\nLevel 3")+
  guides(fill = FALSE)

# Simpler version:
egg::ggarrange(a1, b1, c1, labels = c('A', 'B', 'C'))


# Pearson correlation for each environmental level.

head(df)
hist(log10(df$DR_count))
hist(log10(df$viral_abund))

# All samples
cor.test(log10(df$DR_count), log10(df$viral_abund))

# Run on individual EMPO levels using nest-map-unest workflow 

nested_1 <- df %>% 
  nest(data = -EMPO_level1)

n1<-nested_1 %>% 
  mutate(test = map(data, ~ cor.test(.x$DR_count, .x$viral_abund)),
         tidied = map(test, broom::tidy)) %>%
  unnest(tidied) %>%
  arrange(p.value) %>%
  mutate(EMPO_level = 1) %>%
  rename(environment = EMPO_level1)

# Repeat for EMPO level 2
nested_2 <- df %>% 
  nest(data = -EMPO_level2)

n2<-nested_2 %>% 
  mutate(test = map(data, ~ cor.test(.x$DR_count, .x$viral_abund)),
         tidied = map(test, broom::tidy)) %>%
  unnest(tidied) %>%
  arrange(p.value) %>%
  mutate(EMPO_level = 2) %>%
  rename(environment = EMPO_level2)

# Repeat for EMPO level 3. Check how many obs per group:
keeps<-df %>%
  group_by(EMPO_level3) %>%
  tally() %>%
  filter(n > 5) %>%
  select(EMPO_level3)

keeps<-unique(keeps$EMPO_level3)

nested_3 <- df %>% 
  filter( EMPO_level3 %in% keeps) %>%
  nest(data = -EMPO_level3)

n3<-nested_3 %>% 
  mutate(test = map(data, ~ cor.test(.x$DR_count, .x$viral_abund)),
         tidied = map(test, broom::tidy)) %>%
  unnest(tidied) %>%
  arrange(p.value) %>%
  #mutate( is_sig = ifelse(p.value > 0.05, "no", "yes"))
  mutate(EMPO_level = 3) %>%
  rename(environment = EMPO_level3)

SI_table<-rbind.data.frame(n1, n2, n3)

SI_table<-SI_table %>%
  mutate(adj.p.value = p.adjust(p.value, method = "bonferroni", n = length(p.value))) %>%
  mutate(environment = gsub("\n", " ", environment)) %>%
  select(environment, EMPO_level, estimate, p.value, adj.p.value) %>%
  arrange(EMPO_level, p.value)





