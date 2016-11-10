
library(ggplot2)
library(dplyr)



n <- 50  # number of sequences to do
sl <- 30  # sequence length
bcl <- 4 # bar code length



set.seed(5)


inds <- sample(x = c("Indiv_1", "Indiv_2", "Indiv_3"), size = n, replace = TRUE)
plates <- sample(x = c("Plate_1", "Plate_2", "Plate_3"), size = n, replace = TRUE)


# make a data frame of sequences
df <- data_frame(ID = 1:n) %>%
  mutate(y = sample(ID) * 2,
         xstart = runif(n(), min = 2, max = 120 - sl - 2 * bcl - 2)) %>%
  mutate(bc1end = xstart + bcl,
         bc2start = bc1end + sl,
         bc2end = bc2start + bcl) %>%
  mutate(ind = inds,
         plate = plates)

ggplot(df) +
  geom_segment(mapping = aes(x = xstart, xend = bc1end, y = y, yend = y, colour = ind), size = 1.2) +
  geom_segment(mapping = aes(x = bc1end, xend = bc2start, y = y, yend = y), colour = "darkgrey") +
  geom_segment(mapping = aes(x = bc2start, xend = bc2end, y = y, yend = y, colour = plate), size = 1.2) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) + 
  scale_colour_manual(name = "Barcodes",
                      values = c("red", "orange", "yellow", "violet", "blue", "green"))

ggsave(filename = "figure-creation/output/combi-barcodes.pdf", width = 7, height = 4)
