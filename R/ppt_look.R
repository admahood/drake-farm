# precip history

library(tidyverse)
library(lubridate)

ppt_history <- readxl::read_xlsx("data/Drake Precip All Gages Summary.xlsx",
                                 sheet = 2, skip = 3)[1:12, 1:22] %>%
  mutate(`2006` = as.numeric(`2006`),
         `2015` = as.numeric(`2015`)) %>%
  rename(month = `...1`)

ppt_long <- ppt_history %>%
  pivot_longer(names_to = "year",values_to = "ppt", -month)
