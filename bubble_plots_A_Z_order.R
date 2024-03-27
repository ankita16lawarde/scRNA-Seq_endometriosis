##### redo plots ## A-Z order ##########

## load the saved RData file for saved object to create bubble plot
load("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/scRNA-seq/cell_ranger_count_output/DESeq2/new_cluster_names/merged/cell_cycle_sep_DEG_between_cluster/G2M/S/pathway_bubbleplot/Tuuli_S.RData")

## plot the bubbleplot
ggplot(data = test_df, 
       aes(x = id, y = factor(features.plot, levels = rev(levels(factor(features.plot)))), size = avg.exp)) + 
  #geom_point(aes(colour = pct.exp)) +
  geom_point(data = test_df,  color='grey') +
  geom_point(data = test_df %>% 
               filter(as.numeric(as.character(test_df$padj)) <= 0.05 & test_df$log2FoldChange >= 0.56),
             colour="red") +
  geom_point(data = test_df %>% 
               filter(as.numeric(as.character(test_df$padj)) <= 0.05 & test_df$log2FoldChange <= 0.56),
             colour="blue") +
  theme_bw() + 
  theme(legend.position = "none") +
  ylab("") + 
  xlab("") + 
  ggtitle("Energy Transfer, Membrane Transport")

ggsave("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/scRNA-seq/cell_ranger_count_output/DESeq2/new_cluster_names/merged/cell_cycle_phase_bubbleplots_new_colorcodes_A_Z_order/S_phase/Tuuli.pdf", width = 7.5,
       height = 15)

dev.off()
