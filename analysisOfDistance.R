####################################### Analysis of distance
## Estimating the relations between strains before and after batch correction procedure

############### DAY12

a<-adonis(t(batch_corrected_mat_d12)~ SampleGroup_day12, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ms_data_day12_nonzero_scale) ~ SampleGroup_day12,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day12  21    169158  8055.2  4.9991 0.52235  0.001 ***
#   Residuals          96    154685  1611.3         0.47765           
# Total             117    323844                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Subset data- from each batch
## A subset is used from the scaled data instead of subset the original data set and then scaline
## This is done to preserve the original variation in each datasets as scaling with a subset might have different distributions

## Raw data
d12_batch17_r<-RunDay_day12[RunDay_day12=="17"]
d12_batch17_s<-SampleGroup_day12[names(SampleGroup_day12) %in% names(d12_batch17_r)]#Getting the Sample group and strain names
d12_batch17_data<-ms_data_day12_nonzero_scale[,colnames(ms_data_day12_nonzero_scale) %in% names(d12_batch17_s)]

d12_batch17_data_aod<-adonis(t(d12_batch17_data)~ d12_batch17_s, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(d12_batch17_data) ~ d12_batch17_s, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# d12_batch17_s 11     77775  7070.4   4.719 0.49013  0.001 ***
#   Residuals     54     80907  1498.3         0.50987           
# Total         65    158682                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Batch corrected data
d12_cor_batch17_data<-batch_corrected_mat_d12[,colnames(batch_corrected_mat_d12) %in% names(d12_batch17_s)]

d12_cor_batch17_data_aod<-adonis(t(d12_cor_batch17_data)~ d12_batch17_s, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(d12_cor_batch17_data) ~ d12_batch17_s, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# d12_batch17_s 11     60677  5516.1  4.2697 0.46517  0.001 ***
#   Residuals     54     69762  1291.9         0.53483           
# Total         65    130439                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


############### DAY4

a<-adonis(t(ms_data_day4_nonzero_scale)~ SampleGroup_day4, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ms_data_day4_nonzero_scale) ~ SampleGroup_day4,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day4  21    143357  6826.5   3.605 0.42363  0.001 ***
#   Residuals        103    195041  1893.6         0.57637           
# Total            124    338397                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Subset data- from each batch
## A subset is used from the scaled data instead of subset the original data set and then scaline
## This is done to preserve the original variation in each datasets as scaling with a subset might have different distributions

## Raw data
d4_batch23_r<-RunDay_day4[RunDay_day4=="23"]
d4_batch23_s<-SampleGroup_day4[names(SampleGroup_day4) %in% names(d4_batch23_r)]#Getting the Sample group and strain names
d4_batch23_data<-ms_data_day4_nonzero_scale[,colnames(ms_data_day4_nonzero_scale) %in% names(d4_batch23_s)]

d4_batch23_data_aod<-adonis(t(d4_batch23_data)~ d4_batch23_s, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(d4_batch23_data) ~ d4_batch23_s, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# d4_batch23_s 13     61300  4715.4  2.4771 0.32792  0.001 ***
#   Residuals    66    125636  1903.6         0.67208           
# Total        79    186937                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

d4_cor_batch23_data_nmds <- metaMDS (d4_cor_batch23_data, distance= "euclidean")
d4_batch23_data_nmds <- metaMDS (d4_batch23_data, distance= "euclidean")

## Batch corrected data
d4_cor_batch23_data<-batch_corrected_mat_d4[,colnames(batch_corrected_mat_d4) %in% names(d4_batch23_s)]

d4_cor_batch23_data_aod<-adonis(t(d4_cor_batch23_data)~ d4_batch23_s, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(d4_cor_batch23_data) ~ d4_batch23_s, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# d4_batch23_s 13     56657  4358.2  3.2337 0.3891  0.001 ***
#   Residuals    66     88951  1347.7         0.6109           
# Total        79    145608                 1.0000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pdf("aod_plots.pdf",height=8,width=8)
par(mfrow=c(2,2))
plot(density(d4_batch23_data_aod), xlim=c(0,5), main="DAY4 raw")
plot(density(d12_batch17_data_aod), xlim=c(0,5), main="DAY12 raw")
plot(density(d4_cor_batch23_data_aod), xlim=c(0,5), main="DAY4 corrected")
plot(density(d12_cor_batch17_data_aod), xlim=c(0,5), main="DAY12 corrected")
dev.off()
