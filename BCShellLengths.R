library(readxl);library(dplR);library(ggplot2)

# pull in adjusted, final crossdated shell growth annuli
# sum all growth annuli
# divide individual growth annuli by sum of all growth annuli
# multiply that proportion by the 'shell extent' to get ~width
AxW

# use shell dim to get regression coefficient to translate width to length
#### anchoring regressions 
anchor<-ShellDim %>% 
  group_by(Site.Agg, Species) %>%
  filter(!is.na(Length),!is.na(Width),
         Species %in% c("LORN","LCAR","APLI","QVER"))%>%
  slice(1) %>%
  select(Site.Agg, Species, Length, Width) %>%
  mutate(Length=0,Width=0)

reg<-ShellDim %>% bind_rows(anchor) %>%
  group_by(Site.Agg, Species) %>%
  filter(!is.na(Length),!is.na(Width),
         Species %in% c("LORN","LCAR","APLI","QVER"))%>%
  summarize(Ha=lm(Length~Width)$coefficients[1],
         Hb=lm(Length~Width)$coefficients[2])

# apply regression to build age x length table
AxL<-AxW %>% left_join(reg, by=c("Site.Agg","Species")) %>%
  mutate(est.Length=Ha+est.Width*Hb)

# assess difference between final est.Length and real Length
AxL %>% left_join(ShellDim) %>%
  select(Shell.ID, Site.Agg, Species, Length, est.Length) %>%
  mutate(dif.Length=Length-est.Length)

# vonBertanlanffys
VB

