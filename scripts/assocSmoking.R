comb.res.df %>%
  filter(Batch == "Esteller") %>%
  group_by(Sample_Name, method) %>%
  summarize(n = n()) %>%
  left_join(colData(comb.gset.Esteller) %>% data.frame() %>% select(Sample_Name, msmk),  by = "Sample_Name") %>%
  ggplot(aes(x = msmk, y = n)) + geom_boxplot() +
  theme_bw() + facet_grid(~ method)
  

ind.res.df %>%
  filter(Batch == "Esteller") %>%
  group_by(Sample_Name, method) %>%
  summarize(n = n()) %>%
  left_join(colData(ind.Esteller) %>% data.frame() %>% select(Sample_Name, msmk),  by = "Sample_Name") %>%
  ggplot(aes(x = msmk, y = n)) + geom_boxplot() +
  theme_bw() + facet_grid(~ method)
