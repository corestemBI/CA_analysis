# Corestem project - FACS analysis.proj
# 20220208 Sehwan Chun at Corestem, Inc.
# 2.1. CA efficacy test

#### 1. source Loading ####

load("./data/ca_analysis.image")

#### 2. efficacy test with SARA ####

efficacy_test_with_plot = function(efficacy_file, visit_table, index){
    
    #### calculate SARA slope ####
    
    efficacy_table$V12 = (efficacy_table$V2 - efficacy_table$V1) / ((visit_table$V2 - visit_table$V1) / 28)
    efficacy_table$V23 = (efficacy_table$V3 - efficacy_table$V2) / ((visit_table$V3 - visit_table$V2) / 28)
    efficacy_table$V24 = (efficacy_table$V4 - efficacy_table$V2) / ((visit_table$V4 - visit_table$V2) / 28)
    efficacy_table$V25 = (efficacy_table$V5 - efficacy_table$V2) / ((visit_table$V5 - visit_table$V2) / 28)
    
    #### wilcox test with index #### 
    
    tmp12 = ggpaired(data = efficacy_table, cond1 = "V1", cond2 = "V2", line.color = "gray60", line.size = 1,
                     palette = "aaas", fill = "condition",
                     xlab = "Visit", ylab = paste0(index, " score"),
                     title = paste0("Comparison of ", index, " scores by visit (V1 vs V2)")) +
        stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
        stat_summary(fun.data = function(x) data.frame(y=5,label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
    
    tmp23 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V3", line.color = "gray60", line.size = 1,
                     palette = "aaas",fill = "condition",
                     xlab = "Visit", ylab = paste0(index, " score"),
                     title = paste0("Comparison of ", index, " scores by visit (V2 vs V3)")) +
        stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
        stat_summary(fun.data = function(x) data.frame(y=5,label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
    
    tmp24 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V4", line.color = "gray60", line.size = 1,
                     palette = "aaas",fill = "condition",
                     xlab = "Visit", ylab = paste0(index, " score"),
                     title = paste0("Comparison of ", index, " scores by visit (V2 vs V4)")) +
        stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
        stat_summary(fun.data = function(x) data.frame(y=5,label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
    
    tmp25 = ggpaired(data = efficacy_table[-2,], cond1 = "V2", cond2 = "V5", line.color = "gray60", line.size = 1,
                     palette = "aaas",fill = "condition",
                     xlab = "Visit", ylab = paste0(index, " score"),
                     title = paste0("Comparison of ", index, " scores by visit (V2 vs V5)")) +
        stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
        stat_summary(fun.data = function(x) data.frame(y=5,label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
    
    #### wilcox test with slope ####
    
    tmp1223 = ggpaired(data = efficacy_table, cond1 = "V12", cond2 = "V23", line.color = "gray60", line.size = 1, palette = "aaas", fill = "condition", xlab = "Period", ylab = paste0(index, " slope"), title = paste0("Comparison of ", index, " slope (V1-V2 vs V2-V3)")) +
        stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
        stat_summary(fun.data = function(x) data.frame(y=min((c(efficacy_table$V12,efficacy_table$V23)) - 1),label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
    
    tmp1224 = ggpaired(data = efficacy_table, cond1 = "V12", cond2 = "V24", line.color = "gray60", line.size = 1, palette = "aaas", fill = "condition", xlab = "Period", ylab = paste0(index, " slope"), title = paste0("Comparison of ", index, " slope (V1-V2 vs V2-V4)")) +
        stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
        stat_summary(fun.data = function(x) data.frame(y=min((c(efficacy_table$V12,efficacy_table$V24)) - 1),label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
    
    tmp1225 = ggpaired(data = efficacy_table[-2,], cond1 = "V12", cond2 = "V25", line.color = "gray60", line.size = 1, palette = "aaas", fill = "condition", xlab = "Period", ylab = paste0(index, " slope"), title = paste0("Comparison of ", index, " slope (V1-V2 vs V2-V5)")) +
        stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
        stat_summary(fun.data = function(x) data.frame(y=min(na.omit((c(efficacy_table$V12,efficacy_table$V25))) - 1),label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
    
    #### save plots ####
    
    ggsave(paste0("./output/V12_wilcox_", index, "_paired_boxplot.tiff"), tmp12, width = 8, height = 8, compression = "lzw")
    ggsave(paste0("./output/V23_wilcox_", index, "_paired_boxplot.tiff"), tmp23, width = 8, height = 8, compression = "lzw")
    ggsave(paste0("./output/V24_wilcox_", index, "_paired_boxplot.tiff"), tmp24, width = 8, height = 8, compression = "lzw")
    ggsave(paste0("./output/V25_wilcox_", index, "_paired_boxplot.tiff"), tmp25, width = 8, height = 8, compression = "lzw")
    
    ggsave(paste0("./output/V1223_wilcox_", index, "_paired_boxplot.tiff"), tmp1223, width = 8, height = 8, compression = "lzw")
    ggsave(paste0("./output/V1224_wilcox_", index, "_paired_boxplot.tiff"), tmp1224, width = 8, height = 8, compression = "lzw")
    ggsave(paste0("./output/V1225_wilcox_", index, "_paired_boxplot.tiff"), tmp1225, width = 8, height = 8, compression = "lzw")
}
SARA_efficacy_test = efficacy_test_with_plot(efficacy_file_SARA, visit_table, "SARA")



