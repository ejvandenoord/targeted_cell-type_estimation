

reg_celltype = function(bulk_data,ref_data,boundarycelltype=T) { 
  
  # ref_data = ref_panel;boundarycelltype=F
  plot=F
  library(penalized)
  
  bulk_data       = as.matrix( bulk_data)
  ref_data        = as.matrix( ref_data)
  n_samples       = dim(bulk_data)[1]
  n_sites         = dim(bulk_data)[2]
  n_cellTypes     = dim(ref_data)[2]
  
  props           = matrix(NA,n_samples,n_cellTypes)
  colnames(props) = colnames(ref_data)
  r2              = rep(0,n_samples)
  
  bol_positive    = rep(T,n_cellTypes)
  
  if (boundarycelltype ) print( "Start estimation: Use boundary contraint >= 0") else print("Start estimation: No boundary contraint")
  
  if (plot) {
    #dev.off()
    pdf(paste0(analysis_label[i],"_allCG_",all_CG,"_refpanel.pdf"))
    publ_col=c("#D94D4C","#ECA538","#86AA66","#4CB2D1")
    names(publ_col)= c("red","yellow","green","blue")
    max_y=1
    max_x=max_y}
  
  for(i in 1:n_samples) { # i =1 
    
    if (!boundarycelltype) {
      
      mod    = lm(bulk_data[i,] ~ -1 + ref_data)
      coeffs = coefficients(mod)
      
      if (plot) {
        title =  rownames(bulk_data)[i]
        plot(ref_data[,1],coeffs[1]*bulk_data[i,],xlim=c(0,max_x),ylim=c(0,max_y),pch=16,col=publ_col[1],main=title,ylab="bulk",xlab="ref_data" )
        for (j in 2:n_cellTypes) points(ref_data[,j],coeffs[j]*bulk_data[i,],col=publ_col[j],pch=16)
        
        
        lines(ref_data[,1],coeffs[1]*ref_data[,1],col=publ_col[1])
        for (j in 2:n_cellTypes) lines(ref_data[,j],coeffs[j]*ref_data[,j],col=publ_col[j])
        
        if ( any( grepl("EPIT",rownames(ref_data)))) {
          points(rep(bulk_data[i,2],n_cellTypes),ref_data[2,],col="black",pch=8)
          points(rep(bulk_data[i,8],n_cellTypes),ref_data[8,],col="black",pch=8) }
        
        legend("topleft",legend=colnames(ref_data),pch=rep(16,n_cellTypes),col=publ_col[1:n_cellTypes],bty="n");
        
      }
      
      
      r2[i]        = summary(mod)$r.squared
      props[i,1:n_cellTypes] = coef(mod)[1:n_cellTypes] } else {
        
        sel = complete.cases( cbind(bulk_data[i,],ref_data,row.names = NULL) )
        mod = penalized(bulk_data[i,sel], ~ ref_data[sel,], ~-1,lambda1=0, lambda2=0, positive=bol_positive,trace=F)
        r2[i] = cor( fitted.values(mod), bulk_data[i,sel],use="pairwise.complete.obs")^2
        props[i,1:n_cellTypes] = coef(mod)[1:n_cellTypes]           
        
      }
  }  # for(i in 1:n_samples) 
  
  if (plot) dev.off()
  
  props[is.na(props)] = 0
  nprops              = data.frame(props)
  nprops[nprops < 0]  = 0
  estimates           = nprops/rowSums(nprops) # standardize
  rownames(estimates) = rownames(bulk_data)
  colnames(estimates) = colnames(props)
  
  colnames(nprops)    = paste0(colnames(props),"raw")
  nprops$nprops_max1  = apply(nprops[,colnames(nprops)], 1, max) 
  nprops$nprops_min0  = apply(nprops[,colnames(nprops)], 1, min) 
  nprops$nprops_max1  = round(nprops$nprops_max1 == 1) 
  nprops$nprops_min0  = round(nprops$nprops_min0 == 0) 
  
  results             = list()
  results$estimates   = estimates
  results$qc          = data.frame(nprops,r2,boundarycelltype) 
  
  results
  
} 


create_reference = function(data) {
  
  # data = ref_data
  
  sample_ct = rownames(data)
  sel = grepl('epi',tolower(rownames(data)) )
  sample_ct[!sel] = sapply(strsplit(rownames(data)[!sel],'_'),'[[',2)
  sample_ct[sel]  = 'EPIT'

  celltypes       = unique( sample_ct ) 
  n_samples       = dim(data)[1]
  n_sites         = dim(data)[2]
  n_cellTypes     = length(celltypes)
  
  # nsites by ncell types
  ref_data = matrix(NA,n_sites,n_cellTypes)
  colnames( ref_data ) = celltypes
  rownames( ref_data ) = colnames(data)
  i = 1
  for (i in seq_along(celltypes)) {
    sel_row = sample_ct == celltypes[i]
    ref_data[,i] = colMeans( data[sel_row,] , na.rm = T)
  }
  
   
  clusters             = as.data.frame(cbind(rownames(data),sample_ct) )
  colnames( clusters ) = c("cell_id","clusters")
  rownames( clusters ) = rownames( data )
  clusters             = AnnotatedDataFrame( clusters  )
  
  ref_data_eset        = ExpressionSet(assayData = as.matrix(t(data)), phenoData = clusters)
  
  list(ref_data,ref_data_eset)

}  



duplicate_site_correlations=function(label,methylation_data,outlier_crit) {
  
  # label="25ng"  
  
  sel_samples=grep(label,rownames(methylation_data))
  data=methylation_data[sel_samples,sel_targets]
  sample_counts=table(methylation_data[sel_samples,"sample_ids"])
  n_samples=length(sample_counts)
  
  names(sample_counts)
  
  
  results=matrix(NA,n_targets,7)
  rownames(results)=paste0(label,"_",colnames(methylation_data)[sel_targets])
  colnames(results)=  c("_n_NA","_n_outlier","_mean","_SD","_min","_max","_cor")
  
  j=1
  for (j in 1:n_targets) {
    
    dupl_data=matrix(NA,n_samples,max(sample_counts))
    
    i=1
    for (i in 1:n_samples) dupl_data[i,(1:sample_counts[i])]=data[grep(names(sample_counts)[i],rownames(data)),j]
    n_outlier=sum( is.na(dupl_data) )  
    if (outlier_crit>0) {
      m=mean( dupl_data, na.rm = T)
      s=sd( dupl_data, na.rm = T)
      sel=!(dupl_data > (m+outlier_crit*s)) & !(dupl_data < (m-outlier_crit*s)) 
      k=1
      for (k in 1:length(sel[1,])) dupl_data[!sel[,k],k]=NA
    }  
    
    results[j,1]=sum( is.na(dupl_data) ) 
    results[j,2]= results[j,1]-n_outlier
    results[j,3]=mean( dupl_data, na.rm = T)
    results[j,4]=sd( dupl_data, na.rm = T)
    results[j,5]=min( dupl_data, na.rm = T)
    results[j,6]=max( dupl_data, na.rm = T)
    
    cor=suppressWarnings(cor(dupl_data,use = "pairwise.complete.obs"))
    results[j,7]=mean( cor[lower.tri(cor)],na.rm = T)
    
    
  }
  
  results
}  

parse_sample_names = function(celltype_samples) {
  
  celltypes     = c("CD03","CD3","CD14","CD15","CD19","EPIT","BC","pool")
  samples           = matrix(NA,length(celltype_samples),3)
  colnames(samples) = c("names","type","duplicate")
  rownames(samples) = celltype_samples
  samples[,"names"] = paste0("s",substring(celltype_samples,1,3))
  for (i in seq_along(celltypes)) samples[grepl(celltypes[i],celltype_samples),"type"] = celltypes[i] 
  duplicate_label = c("100D","100_D","_D")  
  for (i in seq_along(duplicate_label)) samples[grepl(duplicate_label[i],celltype_samples),"duplicate"] = "D" 
  samples = samples[order(samples[,"names"],samples[,"type"],samples[,"duplicate"]),]
  samples
}

parse_amplicon_names = function(amplicons,colnames_data,n_short_pos)  {
  
 
  # colnames_data = colnames(data) 
  celltypes     = c("CD03","CD3","CD14","CD15","CD19","EPIT","BC","pool")
  amplicon_type = rep(NA,nrow(amplicons))
  names(amplicon_type) = paste0("chr",amplicons$chr,"_",amplicons$pos )
  for (i in seq_along(celltypes)) amplicon_type[grepl(celltypes[i],amplicons$primersetname)] = celltypes[i]
  amplicon_type
  
#  dim( amplicons )
#  dim( result  )
#  length( short_pos )
  
  
  temp      = strsplit(colnames_data,"_" )
  pos       = unlist(lapply(temp, `[`, 2))
  short_pos = substr(pos,1,n_short_pos)
  result = data.frame( as.integer( gsub("chr","",unlist(lapply(temp, `[`, 1))) ),
                  as.integer(pos), short_pos            )
  rownames(result) = colnames_data
  colnames(result) = c('chr','pos','short_pos')
  target_pos = as.integer( !is.na( match(colnames_data,names(amplicon_type)) ) )
  result      =  data.frame(result,target_pos)
  
  
  temp      = strsplit(names(amplicon_type),"_" )
  chr       = unlist(lapply(temp, `[`, 1))
  pos       = unlist(lapply(temp, `[`, 2))
  short_pos = substr(pos,1,n_short_pos)
  ind = match(paste0("chr",result[,"chr"],"_",result[,"pos"]),names(amplicon_type))
  amplicon_type = amplicon_type[ind]
  
  primersetname = amplicons$primersetname[ind]
  result  =  data.frame(result,amplicon_type,primersetname)
 
  short_pos_string = paste0("chr",result[,"chr"],"_",result[,"short_pos"])
  short_pos_unique = unique( short_pos_string  )
  for ( i in seq_along(short_pos_unique)) { # i = 2

    amplicon_type = result[ short_pos_unique[i] == short_pos_string & result$target_pos==1,'amplicon_type' ][1]
    result[ short_pos_unique[i] == short_pos_string & result$target_pos==0,'amplicon_type' ] = amplicon_type
    
  }
   
  result
  
  
}
