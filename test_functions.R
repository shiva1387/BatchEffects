###############################################
## R functions- batch effect removal ##########
###############################################

#function to compute PCA and perform linear models against runday and strain
compute_pca<-function(dataset,preprocess_method) {
    dataset<-as.data.frame(log(dataset)) #log transform the data using natural log
    if(preprocess_method=="norm") {
      processed_data<-normalize.quantiles(as.matrix(dataset),copy=TRUE)
    }   else  {
      processed_data<-scale(dataset,center=T,scale=T)
      processed_data<-processed_data-min(processed_data)
    }
    pca_results <- princomp(processed_data,cor=F,scores=T) ### IMP: choose quantile normalized or scaled data
    return(pca_results)
  }

#function to compute linear model 
compute_linearModel<-function(results.from.pca,dependent.factor) { #dependent.factor is either RunDay_day4 or 12 (or) SampleGroup_day4 or 12 
    lm_pca_scores<-apply(results.from.pca$loadings,2, function(x) {
      lm_val<-lm(x~ as.factor(dependent.factor)) 
      lm_cor<-summary(lm_val)
      p.val<-anova(lm_val)$'Pr(>F)'[1]
      return(list(lm_cor$r.squared,p.val))
    })
  }

#function to extract r2 value from list containg r2 and p.val returned from linear model
compute.r2.pval<-function(linearmodel_list,r2.pval) {
   if(r2.pval=="r2") { #WARNING:code implictly assumes r2 is in the first column and p.val in the second
     return (sapply(linearmodel_list, function(x){as.numeric(x[1])}))
   } else{
     return (sapply(linearmodel_list, function(x){as.numeric(x[2])}))
   }
  }

#function to compute svd
compute_svd<-function(dataset,preprocess_method,start_pc_comp,end_pc_comp,recursive_pcs=FALSE) {
    dataset<-as.data.frame(log(dataset)) #log transform the data using natural log
    if(preprocess_method=="norm") {
      processed_data<-normalize.quantiles(as.matrix(dataset),copy=TRUE)
    }   else  {
      processed_data<-scale(dataset,center=T,scale=T)
      processed_data<-processed_data-min(processed_data)
    }
    
    if(recursive_pcs) {
      svd_dataset<-svd(processed_data)
      if(end_pc_comp){end_pc<-end_pc_comp} else{end_pc<-ncol(svd_dataset)}
      pc_combi_list<-sapply(start_pc_comp:end_pc,function(x) seq(start_pc_comp:x))
      dataset_rmbatch<-list()
      for(i in 1:length(pc_combi_list))
      {
      svd_dataset$d1<-svd_dataset$d
      svd_dataset$d1[c(unlist(pc_combi_list[i]))]<-0 #Removing the variation caused by runday(id using pca) where the multiple r2 correlation is above 0.5
      dataset_rmbatch1<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v)
      rownames(dataset_rmbatch1)<-rownames(dataset)
      colnames(dataset_rmbatch1)<-names(dataset)
      dataset_rmbatch[[i]]<-dataset_rmbatch1
      }       
    } else{
      dataset_rmbatch<-list()
      svd_dataset<-svd(processed_data)
      svd_dataset$d1<-svd_dataset$d
      svd_dataset$d1[start_pc_comp]<-0 #Removing the variation caused by runday(id using pca) where the multiple r2 correlation is above 0.5
      dataset_rmbatch<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v)
      rownames(dataset_rmbatch)<-rownames(dataset)
      colnames(dataset_rmbatch)<-names(dataset)
    }
    return(dataset_rmbatch)
  }

# compute permutative f test statistics to difentify signigicant features
compute_perm_ftest<-function(dataset,classlabel) {
    classlabel_factor<-as.numeric(as.factor(classlabel))-1
    if(class(dataset) == "list") {
      p.values<-list()
      sig_metab_dataset<-list()
      for(i in 1:length(dataset))
      { data_matrix<-as.data.frame(dataset[i])
        dataset_sig_features<-mt.maxT(data_matrix,classlabel_factor,test="f",side="abs",fixed.seed.sampling="y",B=100000,nonpara="n")
        p.values[[i]]<-dataset_sig_features
        id.sig_dataset <- which(dataset_sig_features$adjp < 0.05 );
        metab_sig<-cbind(data_matrix[id.sig_dataset,],round(dataset_sig_features$adjp[id.sig_dataset],5))
        sig_metab_dataset[[i]]<-metab_sig
      }
     } else { # only a single dataset
      cat ("not a list")
     }
      return(list(p.values,sig_metab_dataset))
  }

#computer list of p.val
compute_pval_list<-function(sigfeat_pval) {
  listofpval<-list()
  for(i in 1:length(sigfeat_pval))
  {
    listofpval[i]<-sigfeat_pval[[i]][4]
  }
  sigfeat_pval_pvaldf<-do.call(cbind.data.frame, listofpval)
}

#compute no of sig features
compute_no_features<-function(sig_feat) {
  no_of_sig_features<-rep(NA,length(sig_feat))
  for(i in 1:length(sig_feat))
  {
    no_of_sig_features[i]<-nrow(sig_feat[[i]][1])
  }
  return(no_of_sig_features)
}
  
###############################################
#day4

fit_day4_mzbysam_princomp<-compute_pca(ms_data_day4_nonzero,"scale")
lm_pca_scores_strain_day4_nonzero_loadings<-compute_linearModel(fit_day4_mzbysam_princomp,SampleGroup_day4)
lm_pca_scores_strain_day4_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_strain_day4_nonzero_loadings,"r2")
lm_pca_scores_strain_day4_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_strain_day4_nonzero_loadings,"pval")

lm_pca_scores_runday4_nonzero_loadings<-compute_linearModel(fit_day4_mzbysam_princomp,RunDay_day4)
lm_pca_scores_runday4_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_runday4_nonzero_loadings,"r2")
lm_pca_scores_runday4_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_runday4_nonzero_loadings,"pval")

svd_day4_nonzero <- compute_svd(ms_data_day4_nonzero,"scale",1,121,TRUE)
#str(svd_day4_nonzero)
day4_nonzero_sigfeat_r<-compute_perm_ftest(svd_day4_nonzero,RunDay_day4)
day4_nonzero_sigfeat_r_pval<-day4_nonzero_sigfeat_r[[1]]
day4_nonzero_sigfeat_r_matrix<-day4_nonzero_sigfeat_r[[2]]
day4_nonzero_sigfeat_r_pvaldf<-compute_pval_list(day4_nonzero_sigfeat_r_pval)
day4_nonzero_sigfeat<-compute_no_features(day4_nonzero_sigfeat_r_matrix)

###################
#### day12
###################

fit_day12_mzbysam_princomp<-compute_pca(ms_data_day12_nonzero,"scale")

lm_pca_scores_strain_day12_nonzero_loadings<-compute_linearModel(fit_day12_mzbysam_princomp,SampleGroup_day12)
lm_pca_scores_strain_day12_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_strain_day12_nonzero_loadings,"r2")
lm_pca_scores_strain_day12_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_strain_day12_nonzero_loadings,"pval")

lm_pca_scores_runday12_nonzero_loadings<-compute_linearModel(fit_day12_mzbysam_princomp,RunDay_day12)
lm_pca_scores_runday12_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_runday12_nonzero_loadings,"r2")
lm_pca_scores_runday12_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_runday12_nonzero_loadings,"pval")

svd_day12_nonzero <- compute_svd(ms_data_day12_nonzero,"scale",1,124,TRUE)
#str(svd_day12_nonzero)
day12_nonzero_sigfeat_r<-compute_perm_ftest(svd_day12_nonzero,RunDay_day12)
day12_nonzero_sigfeat_r_pval<-day12_nonzero_sigfeat_r[[1]]
day12_nonzero_sigfeat_r_matrix<-day12_nonzero_sigfeat_r[[2]]
day12_nonzero_sigfeat_r_pvaldf<-compute_pval_list(day12_nonzero_sigfeat_r_pval)
day12_nonzero_sigfeat<-compute_no_features(day12_nonzero_sigfeat_r_matrix)

matplot(day4_nonzero_sigfeat_r_pvaldf,type='l')
matplot(day12_nonzero_sigfeat_r_pvaldf,type='l')
plot(1:121,day12_nonzero_sigfeat,type='l',col='blue')
plot(1:121,day4_nonzero_sigfeat,type='l',col='red')


