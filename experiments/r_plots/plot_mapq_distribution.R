library("ggplot2")
# df_mapq_dist <- read.delim("21-mapq.txt", header=F)
df_mapq_dist <- read.delim("/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow/wg_hg38_2/experiments/NA12878/wg-GRC.sam.mapq", header=F)
head(df_mapq_dist)

p_hist_mapq <- ggplot(df_mapq_dist, aes(V1)) +
  scale_y_continuous(trans = 'log10', limits=c(1,NA)) + theme_bw() + geom_histogram(bins=43, closed="left") +
  xlab("MAPQ") + ylab("Count")
ggsave("fig_hist_count_mapq.pdf", p_hist_mapq, width = 12, height = 8)

p_freq_mapq <- ggplot(df_mapq_dist, aes(V1, stat(density))) +
  theme_bw() + geom_freqpoly(bins=43, closed="left") +
  xlab("MAPQ") + ylab("Density")
ggsave("fig_hist_freq_mapq.pdf", p_freq_mapq, width = 12, height = 8)

p_cumcount_mapq <- ggplot(df_mapq_dist, aes(V1)) +
  stat_ecdf(geom="step") +
  theme_bw() +
  xlab("MAPQ") + ylab("Cumulative frequency")
ggsave("fig_cumcount_mapq.pdf", p_cumcount_mapq, width = 12, height = 8)
