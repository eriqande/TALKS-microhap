
library(ggplot2)
library(dplyr)



n <- 150  # number of sequences to do
sl <- 20  # sequence length
bcl <- 3 # bar code length



set.seed(5)


inds <- sample(x = c("Indiv_1", "Indiv_2", "Indiv_3"), size = n, replace = TRUE)
plates <- sample(x = c("Plate_1", "Plate_2", "Plate_3"), size = n, replace = TRUE)
amplicons <- sample(c("Amplicon_1", "Amplicon_2"), size = n, replace = TRUE)


# make a data frame of sequences
df <- data_frame(ID = 1:n) %>%
  mutate(y = sample(ID) * 2,
         xstart = runif(n(), min = 2, max = 120 - sl - 2 * bcl - 2)) %>%
  mutate(bc1end = xstart + bcl,
         bc2start = bc1end + sl,
         bc2end = bc2start + bcl) %>%
  mutate(ind = inds,
         plate = plates,
         amplicon = amplicons)


ggplot(df) +
  geom_segment(mapping = aes(x = xstart, xend = bc1end, y = y, yend = y, colour = ind), size = 1.0) +
  geom_segment(mapping = aes(x = bc1end, xend = bc2start, y = y, yend = y, colour = amplicon)) +
  geom_segment(mapping = aes(x = bc2start, xend = bc2end, y = y, yend = y, colour = plate), size = 1.0) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) + 
  scale_colour_manual(name = "Amplicons\n and\nBarcodes",
                      values = c("lightgray", "black", "red", "orange", "yellow", "violet", "blue", "green"))




ggsave(filename = "figure-creation/output/gtseq-soup.pdf", width = 7, height = 4)



# then pull them all together into two different amplicons  
ampdb <- data_frame(
  label = c("Reference 1", "Reference 2"),
  y = 0,
  rectxstart = c(20 + bcl, 75 + bcl)
) %>%
  mutate(rectxend = rectxstart + sl,
         textx = rectxstart + 0.5,
         recty = -2)

df2 <- df %>%
  mutate(xstart = ifelse(amplicon == "Amplicon_1", 20, 75)) %>%
  mutate(bc1end = xstart + bcl,
         bc2start = bc1end + sl,
         bc2end = bc2start + bcl)


ggplot(df2) +
  geom_segment(mapping = aes(x = xstart, xend = bc1end, y = y, yend = y, colour = ind), size = 1.0) +
  geom_segment(mapping = aes(x = bc1end, xend = bc2start, y = y, yend = y, colour = amplicon)) +
  geom_segment(mapping = aes(x = bc2start, xend = bc2end, y = y, yend = y, colour = plate), size = 1.0) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) + 
  scale_colour_manual(name = "Amplicons\n and\nBarcodes",
                      values = c("lightgray", "black", "red", "orange", "yellow", "violet", "blue", "green")) +
  xlim(0,120) +
  geom_rect(data = ampdb, mapping = aes(xmin = rectxstart, xmax = rectxend, ymax = recty, ymin = recty - 18)) +
  geom_text(data = ampdb, mapping = aes(x = textx, y = recty - 9, label = label), colour = "white", hjust = 0, size = 8)
  



ggsave(filename = "figure-creation/output/gtseq-amps.pdf", width = 14, height = 8)



# then sort them by barcode and plot them together. This is sort of klugie, but I basically want each 
# individual to take up a height of 16.  
indht <- 16
ygroups <- df2 %>%
  group_by(ind, plate) %>%
  tally() %>%
  ungroup() %>%
  mutate(group_y = indht * 0:(n() - 1)) %>%
  select(-n)

df3 <- df2 %>%
  left_join(., ygroups) %>%
  arrange(amplicon, ind, plate) %>%
  group_by(amplicon, ind, plate) %>%
  mutate(y = group_y + 1:n())



g <- ggplot(df3) +
  geom_segment(mapping = aes(x = xstart, xend = bc1end, y = y, yend = y, colour = ind), size = 1.0) +
  geom_segment(mapping = aes(x = bc1end, xend = bc2start, y = y, yend = y, colour = amplicon)) +
  geom_segment(mapping = aes(x = bc2start, xend = bc2end, y = y, yend = y, colour = plate), size = 1.0) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) + 
  scale_colour_manual(name = "Amplicons\n and\nBarcodes",
                      values = c("lightgray", "black", "red", "orange", "yellow", "violet", "blue", "green")) +
  xlim(0,120) +
  geom_hline(yintercept = 0:8 * indht, linetype = "dashed")  +
  geom_rect(data = ampdb, mapping = aes(xmin = rectxstart, xmax = rectxend, ymax = recty, ymin = recty - 9)) +
  geom_text(data = ampdb, mapping = aes(x = textx, y = recty - 4.5, label = label), colour = "white", hjust = 0, size = 8.0)



ggsave(filename = "figure-creation/output/gtseq-demultiplexed.pdf", width = 14, height = 8)


# now, we are going to want to create haplotypes for these individuals to carry.
haps <- read.table("figure-creation/inputs/haplos.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  tbl_df()

# and now we are going to want to assign haplotypes to each individual randomly
assign_haps <- function(amp) {
  if(amp[1] == "Amplicon_1") {
    nHap <- 5
  } else {
    nHap <- 2
  }
  n <- length(amp)
  
  nh1 <- rbinom(1, size = n, prob = 0.5)
  haps <- sort(c(sample(1:nHap, size = 1), sample(1:nHap, size = 1)))
  rep(haps, c(nh1, n - nh1)) 
}

df4 <- df3  %>%
  group_by(amplicon, ind, plate) %>%
  mutate(haplotype = assign_haps(amplicon)) %>%
  ungroup() %>%
  left_join(haps)  %>%
  mutate(snp_x = bc1end + pos)

g + 
  geom_point(data = df4, mapping = aes(x = snp_x, y = y), colour = "red", size = 0.85)

ggsave(filename = "figure-creation/output/gtseq-snps.pdf", width = 14, height = 8)
