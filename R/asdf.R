source("R/drake_data_prep.R")

ggplot(xdata, aes(x = soil_moisture_pre_ma_30cm, fill=strip_type)) +
  geom_histogram()

ggplot(xdata, aes(x = soil_moisture_pre_son_30cm, fill=strip_type)) +
  geom_histogram()

ggplot(soil_moisture_summaries, aes(x=mean_soil_moisture)) +
  geom_histogram() +
  facet_wrap(~depth)
