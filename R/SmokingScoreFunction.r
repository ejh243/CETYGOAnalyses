smokingscore<-function(betas){
  
	load("RData/SmokingScoreRefData.rda")
	#contains Illig_data, Illig_data_up and Illig_data_down
  
	#subsetting own data 
	#SABRE_data_down
	select <- rownames(betas) %in% Illig_data_down$cpgs
	A_down <- subset(betas, select =="TRUE")
  
	#SABRE_data_up
	select <- rownames(betas) %in% Illig_data_up$cpgs
	A_up <- subset(betas, select =="TRUE")
  
	#sort SABRE data by Cpg name
	A_up <- A_up[order(rownames(A_up)),]
	A_down <- A_down[order(rownames(A_down)),]
		
	#match Illig data by by Cpg name
	Illig_data_up<-Illig_data_up[match(rownames(A_up), Illig_data_up$cpgs),]
	Illig_data_down<-Illig_data_down[match(rownames(A_down), Illig_data_down$cpgs),]
  
	#as some outliers have been removed and replaced with NAs need to handle missing values
	matrix_up_A<- matrix(nrow=nrow(Illig_data_up), ncol=ncol(A_up))
	for (i in 1:ncol(A_up)){
	  matrix_up_A[,i]<- (A_up[,i])-(Illig_data_up$reference_never_median_beta_all)}
	colnames(matrix_up_A)<- colnames(A_up)
	rownames(matrix_up_A)<- Illig_data_up$cpgs

	#Calculate scores - ##UP##
	#calculate scores for each individual in the dataset
	scores_up_A<- as.numeric(rep(NA,ncol(A_up)))
	for (i in 1:ncol(A_up)){
	  scores_up_A[i]<-sum(matrix_up_A[,i]*Illig_data_up$weights)}

	  
	#Calculate diffs between SABRE beta values and the reference for each site - ##DOWN###
	matrix_down_A<- matrix(nrow=nrow(Illig_data_down), ncol=ncol(A_down))
	for (i in 1:ncol(A_down)){
	  matrix_down_A[,i]<- (Illig_data_down$reference_never_median_beta_all)-(A_down[,i])}
	colnames(matrix_down_A)<- colnames(A_down)
	rownames(matrix_down_A)<- Illig_data_down$cpgs

	#Calculate scores - ##DOWN##
	#calculate scores for each individual in the dataset
	scores_down_A<- as.numeric(rep(NA,ncol(A_down)))
	for (i in 1:ncol(A_down)){
	  scores_down_A[i]<-sum(matrix_down_A[,i]*Illig_data_down$weights)}

	##combine scores
	scores_combined_A <- scores_up_A + scores_down_A
  
  return(scores_combined_A)
}