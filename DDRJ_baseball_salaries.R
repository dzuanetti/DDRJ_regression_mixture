#
library(mvtnorm)
library(MCMCpack)
library(compiler)
library(stringr)
#
####### Necessary functions
#
rDiscreta <- function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  return(val)}
#
rDiric = function(a){
  X<-rgamma(length(a),a,1)
  Y<-X/sum(X)
  return(Y)}
#
posterior_sigma = function(prior_a, prior_b, res){
  alpha <- (length(res)/2) + prior_a
  beta <- (sum(res^2)/2) + prior_b
  sigma2 <- 1/(rgamma(1,alpha,beta))
  dens<-dgamma(1/sigma2,alpha,beta,log=TRUE)
  return(c(sigma2,dens))}
#
log_posterior_sigma = function(prior_a, prior_b, res, sigma2_current){
  alpha<-(length(res)/2) + prior_a
  beta<-(sum(res^2)/2) + prior_b
  sigma2_post<-dgamma(1/sigma2_current,alpha,beta,log=TRUE)
  return(sigma2_post)}
#
posterior_beta_mult = function(current_sigma2, Y, X, sigma_beta, mu_beta){
	if (length(c(Y))>=1){
		sigma_beta_inv<-solve(sigma_beta)
		Xtrans<-t(X)
		XlX<-Xtrans%*%X
		var_beta<-solve(XlX+sigma_beta_inv*current_sigma2)
		media_beta<-var_beta%*%(Xtrans%*%Y+((sigma_beta_inv*current_sigma2)%*%mu_beta))
		beta_vig<-rmvnorm(1,mean=c(media_beta), sigma=(var_beta*current_sigma2), pre0.9_9994 = TRUE)
		beta_post<-dmvnorm(c(beta_vig), mean = media_beta, sigma = (var_beta*current_sigma2), log = TRUE)} else{
		beta_vig<-rmvnorm(1,mean=c(mu_beta), sigma=sigma_beta, pre0.9_9994 = TRUE)
		beta_post<-dmvnorm(c(beta_vig), mean = mu_beta, sigma = sigma_beta, log = TRUE)}
	return(list(beta_vig,beta_post))}
#
log_posterior_beta_mult = function(current_sigma2, Y, X, current_beta, sigma_beta, mu_beta){
	if (length(Y)>=1){
		sigma_beta_inv<-solve(sigma_beta)
		Xtrans<-t(X)
		XlX<-Xtrans%*%X
		var_beta<-solve(XlX+sigma_beta_inv*current_sigma2)
		media_beta<-var_beta%*%(Xtrans%*%Y+((sigma_beta_inv*current_sigma2)%*%mu_beta))
		beta_post<-dmvnorm(c(current_beta), mean = media_beta, sigma = (var_beta*current_sigma2), log = TRUE)} else beta_post<-dmvnorm(c(current_beta), mean = mu_beta, sigma = sigma_beta, log = TRUE)
	return(beta_post)}
#
latent_generator = function(tam,n_grupos){
  w_aux<-1/n_grupos
  v_aux<-rep(w_aux,n_grupos)
  s<-NULL
  for(i in 1:tam){
    s[i]<-rDiscreta(v_aux)}
  return(s)}
#
latent_updater<-function(y, x, w, beta1, sigma21, K){
  mu<-NULL
  aux<-NULL
  s<-NULL
  for(i in 1:length(y)){
    for(k in 1:K){
      mu[k]<-x[i,]%*%beta1[,k]
      aux[k]<-dnorm(y[i],mean=mu[k],sd=sqrt(sigma21[k]))*w[k]}
    prob<-aux/sum(aux)
    s[i]<-rDiscreta(prob)}
  return(s)}
#
prob_split_merge = function(K, Kmax){
  p <- NULL
  if(K==1){ 
    p <- c(1,0)
  } else if(K==Kmax){
    p <- c(0,1)  
  } else {p <- c(.5,.5)}
  return(p)}
#
log_likelihood_function<-function(pesos, n_k, sigma2, Sla, res){
  tamK<-length(pesos)
  aux <- NULL
  for (i in 1:tamK) aux[i]<-sum(res[Sla==i]^2)/(2*sigma2[i])
  r <- sum(n_k*log(pesos))-sum(log(2*pi*sigma2)*(n_k/2))-sum(aux)
  return(r)}
#
choose_components_merge<-function(K){
  comp_vector_merge <- sample(1:K,2,replace=FALSE)  
  return(sort(comp_vector_merge))}
#
####### reading the baseball salaries data set - you can analyze other data sets
#
# Here, we need to change for a folder in you computer where the data file is available
dados<-read.table(file="/Users/Daiane/Downloads/baseball.csv",sep=',',header=TRUE)
for (i in 2:13) dados[,i]<-(dados[,i]-mean(dados[,i]))/sd(dados[,i])
Y<-log(dados[,1])
X<-dados[,-1]
n<-nrow(X)
#
# including the interaction covariates in the data set
#
X<-cbind(X,X[,1]*X[,13],X[,1]*X[,14],X[,1]*X[,15],X[,1]*X[,16],X[,3]*X[,13],X[,3]*X[,14],X[,3]*X[,15],X[,3]*X[,16],X[,7]*X[,13],X[,7]*X[,14],X[,7]*X[,15],X[,7]*X[,16],X[,8]*X[,13],X[,8]*X[,14],X[,8]*X[,15],X[,8]*X[,16])
X<-round(as.matrix(cbind(rep(1,n),X)),4)
dim(X)
#
caminho<-"/Users/Daiane/Downloads/" # Again, choose another folder available in your computer where the intermediate files of the MCMC procedure will be recorded
tam<-length(Y)
burnin<-100 # choose the most suitable number for your run
Nfinal<-100 # choose the most suitable number for your run
saltos<-1 # choose the most suitable number for your run
Kmax<-30
Ntotal<-burnin+Nfinal*saltos
pcov<-ncol(X)
#
# hyperparameters specification
shape_prior_var<-0.1
rate_prior_var<-0.1
mean_prior_beta<-0
var_prior_beta<-100
gamma_w<-1
#
#### starting point specification
Kini<-5
n_k<-NULL
res<-NULL
sigma2<-NULL  
K<-Kini
set.seed(100)
S<-latent_generator(length(Y),K)
for (k in 1:K) n_k[k]<-sum(S==k)
betas<-matrix(0,nrow=pcov,ncol=K,byrow=TRUE)
for(i in 1:tam) res[i]<-Y[i]-(X[i,]%*%betas[,S[i]])
for(k in (1:K)) sigma2[k]<-posterior_sigma(shape_prior_var,rate_prior_var,res[S==k])[1]
#
#### starting DDRJ iterations - main block
#
enableJIT(3)
for (iter in (1:Ntotal)){
	cat("iteration = ", iter,"K = ", K, "\n")
	# A) Updating the parameters of the current model via Gibbs Sampling
	for (k in 1:K){
		if (n_k[k]>0){
			betas[,k]<-posterior_beta_mult(sigma2[k],matrix(Y[S==k],ncol=1),matrix(X[S==k,],ncol=pcov),diag(var_prior_beta,ncol=pcov,nrow=pcov),matrix(rep(mean_prior_beta,pcov),ncol=1))[[1]]
 			res[S==k]<-Y[S==k]-(X[S==k,]%*%betas[,k])} else {
			betas[,k]<-posterior_beta_mult(sigma2[k],Y[S==k],X[S==k,],diag(var_prior_beta,ncol=pcov,nrow=pcov),matrix(rep(mean_prior_beta,pcov),ncol=1))[[1]]}
		sigma2[k]<-posterior_sigma(shape_prior_var, rate_prior_var, res[S==k])[1]}
	w<-rDiric(n_k+gamma_w)
	S<-latent_updater(Y,X,w,betas,sigma2,K)
	for (k in 1:K) n_k[k]<-sum(S==k)
	#
	#############################  
	# Updating K    
	#########################################################
	#### Step 1) Choosing for a split or a merge movemente
	#########################################################
	Psm<-prob_split_merge(K, Kmax)
	Dsm<-rDiscreta(Psm) # return 1 if split or 2 if merge
	if(Dsm == 1){ # steps in a split movement
		Ps<-rep(1/K,K)
		Dcomp<-rDiscreta(Ps)
		S_new<-S
		for (i in 1:tam){
			if (S_new[i]==Dcomp){
				auxi<-runif(1)
				if (auxi>0.5) S_new[i]<-K+1}}
		K_new<-K+1     
		n_k_new<-NULL
		for (k in 1:K_new) n_k_new[k]<-sum(S_new==k)
		w_new<-rDiric(n_k_new+gamma_w)
		e_new<-res
		#
		## sampling parameters for new components
		#
		betas_new<-cbind(betas, rep(0, length(betas[,1]))) 
		betas_inf1<-posterior_beta_mult(sigma2[Dcomp],matrix(Y[S_new==Dcomp],ncol=1),matrix(X[S_new==Dcomp,],ncol=pcov),diag(var_prior_beta,ncol=pcov,nrow=pcov),matrix(rep(mean_prior_beta,pcov),ncol=1))
		betas_inf2<-posterior_beta_mult(sigma2[Dcomp],matrix(Y[S_new==K_new],ncol=1),matrix(X[S_new==K_new,],ncol=pcov),diag(var_prior_beta,ncol=pcov,nrow=pcov),matrix(rep(mean_prior_beta,pcov),ncol=1))
		betas_new[,Dcomp]<-betas_inf1[[1]]
		betas_new[,K_new]<-betas_inf2[[1]]
		e_new[S_new==Dcomp]<-Y[S_new==Dcomp]-(X[S_new==Dcomp,]%*%betas_new[,Dcomp])
		e_new[S_new==K_new]<-Y[S_new==K_new]-(X[S_new==K_new,]%*%betas_new[,K_new])
		sigma2_new<-c(sigma2,0) 
		sigma2_inf1<-posterior_sigma(shape_prior_var, rate_prior_var, e_new[S_new==Dcomp])
		sigma2_inf2<-posterior_sigma(shape_prior_var, rate_prior_var, e_new[S_new==K_new])	 
		sigma2_new[Dcomp]<-sigma2_inf1[1]
		sigma2_new[K_new]<-sigma2_inf2[1]
		#
		## calculating the acceptance probability of split
		#
		P1_num<-log_likelihood_function(w_new, n_k_new, sigma2_new, S_new, e_new)
		P1_denom<-log_likelihood_function(w, n_k, sigma2, S, res)
		# 
		P2_num<-log(ddirichlet(w_new, rep(1,length(w_new))))+
	   dgamma(1/sigma2_new[Dcomp], shape_prior_var, rate_prior_var,log=TRUE)+
       dgamma(1/sigma2_new[K_new], shape_prior_var, rate_prior_var,log=TRUE)+
       dmvnorm(c(betas_new[,Dcomp]), mean = matrix(rep(mean_prior_beta,pcov),ncol=1), sigma = diag(var_prior_beta,ncol=pcov,nrow=pcov), log = TRUE)+
       dmvnorm(c(betas_new[,K_new]), mean = matrix(rep(mean_prior_beta,pcov),ncol=1), sigma = diag(var_prior_beta,ncol=pcov,nrow=pcov), log = TRUE)
		P2_denom<-log(ddirichlet(w, rep(1,length(w))))+
       dgamma(1/sigma2[Dcomp], shape_prior_var, rate_prior_var,log=TRUE)+
       dmvnorm(c(betas[,Dcomp]), mean = matrix(rep(mean_prior_beta,pcov),ncol=1), sigma = diag(var_prior_beta,ncol=pcov,nrow=pcov), log = TRUE)
		#        
		P3_num<-log(prob_split_merge(K_new,Kmax)[2])+
       log(2/(K_new*(K_new-1)))+
       log(ddirichlet(w,n_k+gamma_w))+
       log_posterior_sigma(shape_prior_var, rate_prior_var, e_new[S==Dcomp], sigma2[Dcomp])+
	   log_posterior_beta_mult(sigma2[Dcomp], matrix(Y[S==Dcomp],ncol=1), matrix(X[S==Dcomp,],ncol=pcov), betas[,Dcomp], diag(var_prior_beta,ncol=pcov,nrow=pcov), matrix(rep(mean_prior_beta,pcov),ncol=1))
		P3_denom<-log(Psm[1])+log(Ps[Dcomp])+
       n_k[Dcomp]*log(1/2)+
       log(ddirichlet(w_new, n_k_new+gamma_w))+
       betas_inf1[[2]]+betas_inf2[[2]]+
       sigma2_inf1[2]+sigma2_inf2[2]
		prob_ace<-exp((P1_num+P2_num+P3_num)-(P1_denom+P2_denom+P3_denom))
		} else {# steps in a merge movement
		comp_merge<-choose_components_merge(K)
		S_new<-S
		S_new[S==comp_merge[2]]<-comp_merge[1]
		S_new[S>comp_merge[2]]<-S[S>comp_merge[2]]-1
		K_new<-K-1
		n_k_new<-NULL
		for (k in 1:K_new) n_k_new[k]<-sum(S_new==k)
		betas_new<-matrix(betas[,-comp_merge[2]],ncol=K_new,nrow=pcov)      
		sigma2_new<-sigma2[-comp_merge[2]]
		#
		## sampling parameters for new components
		#
		w_new<-rDiric(n_k_new+gamma_w)
		sigma2_inf<-posterior_sigma(shape_prior_var, rate_prior_var, res[S_new==comp_merge[1]])
		sigma2_new[comp_merge[1]]<-sigma2_inf[1]
		betas_inf<-posterior_beta_mult(sigma2_new[comp_merge[1]], matrix(Y[S_new==comp_merge[1]],ncol=1),  matrix(X[S_new==comp_merge[1],],ncol=pcov), diag(var_prior_beta,ncol=pcov,nrow=pcov), matrix(rep(mean_prior_beta,pcov),ncol=1))
		betas_new[,comp_merge[1]]<-betas_inf[[1]]
		e_new<-res
		e_new[S_new==comp_merge[1]]<-Y[S_new==comp_merge[1]]-(X[S_new==comp_merge[1],]%*%betas_new[,comp_merge[1]])
		#
		## calculating the acceptance probability of split
		#
		P1_num<-log_likelihood_function(w_new, n_k_new, sigma2_new, S_new, e_new)
		P1_denom<-log_likelihood_function(w, n_k, sigma2, S, res)
		#
		P2_num<-log(ddirichlet(w_new, rep(1,length(w_new))))+
       dgamma(1/sigma2_new[comp_merge[1]], shape_prior_var, rate_prior_var,log=TRUE)+
       dmvnorm(c(betas_new[,comp_merge[1]]), mean = matrix(rep(mean_prior_beta,pcov),ncol=1), sigma = diag(var_prior_beta,ncol=pcov,nrow=pcov), log = TRUE)
		P2_denom<-log(ddirichlet(w, rep(1,length(w))))+
       dgamma(1/sigma2[comp_merge[1]], shape_prior_var, rate_prior_var,log=TRUE)+
       dgamma(1/sigma2[comp_merge[2]], shape_prior_var, rate_prior_var,log=TRUE)+
       dmvnorm(c(betas[,comp_merge[1]]), mean = matrix(rep(mean_prior_beta,pcov),ncol=1), sigma = diag(var_prior_beta,ncol=pcov,nrow=pcov), log = TRUE)+
       dmvnorm(c(betas[,comp_merge[2]]), mean = matrix(rep(mean_prior_beta,pcov),ncol=1), sigma = diag(var_prior_beta,ncol=pcov,nrow=pcov), log = TRUE)
       #
       P3_num<-log(prob_split_merge(K_new, Kmax)[1])+
       log(1/K_new)+
       n_k_new[comp_merge[1]]*log(1/2)+
       log(ddirichlet(w, n_k+gamma_w))+
       log_posterior_beta_mult(sigma2_new[comp_merge[1]], matrix(Y[S==comp_merge[1]],ncol=1), matrix(X[S==comp_merge[1],],ncol=pcov), betas[,comp_merge[1]], diag(var_prior_beta,ncol=pcov,nrow=pcov), matrix(rep(mean_prior_beta,pcov),ncol=1))+
       log_posterior_beta_mult(sigma2_new[comp_merge[1]], matrix(Y[S==comp_merge[2]],ncol=1), matrix(X[S==comp_merge[2],],ncol=pcov), betas[,comp_merge[2]], diag(var_prior_beta,ncol=pcov,nrow=pcov), matrix(rep(mean_prior_beta,pcov),ncol=1))+
	   log_posterior_sigma(shape_prior_var, rate_prior_var, res[S==comp_merge[1]], sigma2[comp_merge[1]])+
       log_posterior_sigma(shape_prior_var, rate_prior_var, res[S==comp_merge[2]], sigma2[comp_merge[2]])
		P3_denom<-log(Psm[2])+
       log(2/(K*(K-1)))+
       log(ddirichlet(w_new, n_k_new+gamma_w))+
       sigma2_inf[2]+
	   betas_inf[[2]]
	   prob_ace<-exp((P1_num+P2_num+P3_num)-(P1_denom+P2_denom+P3_denom))}
	#
	## evaluating the acceptance of the candidate move
	#
	log_vero<-P1_denom
	aux2<-runif(1)
	if (aux2<prob_ace){
		S<-S_new 
		sigma2<-sigma2_new 
		n_k<-n_k_new
		betas<-betas_new 
		res<-e_new
		K<-K_new
		w<-w_new
		log_vero<-P1_num}
	#
	## recording the MCMC in intermediate files
	#
   if (iter>burnin & iter%%saltos==0){     
     cat('',round(prob_ace,3),file=paste(caminho,"prob_ace.txt",sep=""),append=T)
     cat('',K,file=paste(caminho,"K.txt",sep=""),append=T)
     cat('',round(betas,3),file=paste(caminho,"Betas.txt",sep=""),append=T)
     cat('',round(sigma2,3),file=paste(caminho,"Sigmas2.txt",sep=""),append=T)
     cat('',round(w,3),file=paste(caminho,"Pesos.txt",sep=""),append=T)
     cat('',S,file=paste(caminho,"S.txt",sep=""),append=T)
     cat('',log_vero,file=paste(caminho,"log_vero.txt",sep=""),append=T)}
}
#
#
########
### analyzing the results
########
iter<-Nfinal
pcov<-ncol(X)
#
######
# convergence
######
log_vero<-scan(file=paste(caminho,"log_vero.txt",sep=""))
plot(log_vero,type='l')
library(coda)
log_v<-mcmc(log_vero)
effectiveSize(log_v)
geweke.diag(log_v)
autocorr.plot(log_v, lag.max=20)
#
######
# results
######
K<-scan(file=paste(caminho,"K.txt",sep=""))
Sj<-scan(file=paste(caminho,"S.txt",sep=""))
variancia<-scan(file=paste(caminho,"Sigmas2.txt",sep=""))
baleat<-scan(file=paste(caminho,"Betas.txt",sep=""))
w<-scan(file=paste(caminho,"Pesos.txt",sep=""))
#
plot(K)
round(table(K)/iter*100,1)
#
## rearranging results into matrices
#
m<-length(Y)
S<-matrix(Sj,ncol=m,nrow=iter,byrow=TRUE)
b_aleat<-matrix(0,nrow=iter,ncol=(max(K)*pcov))
var_m<-matrix(0,nrow=iter,ncol=max(K))
w_m<-matrix(0,nrow=iter,ncol=max(K))
obs<-1
obs2<-1
obs3<-1
for (i in 1:iter){
  if (K[i]>0){
    b_aleat[i,1:(K[i]*pcov)]<-baleat[obs2:(obs2+(K[i]*pcov)-1)]
    var_m[i,1:K[i]]<-variancia[obs:(obs+K[i]-1)]
    w_m[i,1:K[i]]<-w[obs:(obs+K[i]-1)]}
  obs<-obs+K[i]
  obs2<-obs2+(K[i]*pcov)}
#
head(w_m)
head(var_m)
head(b_aleat)
#
#### calculating the co-clustering probability
Sj.j<-S
prob.eq<-matrix(0,nrow=ncol(Sj.j),ncol=ncol(Sj.j))
for (i in 1:ncol(Sj.j)){
	for (j in 1:ncol(Sj.j)){
		prob.eq[i,j]<-round(sum(Sj.j[,i]==Sj.j[,j])/length(Sj.j[,i]),5)*100}}
#
# defining the final groups according with the chosen cuttof
#
thresh<-0.50*100 
clust_f<-c(1,rep(0,(ncol(Sj.j)-1)))
for (i in 2:ncol(Sj.j)){
	if (max(prob.eq[i,1:(i-1)])>thresh) clust_f[i]<-clust_f[which(prob.eq[i,1:(i-1)]==max(prob.eq[i,1:(i-1)]))[1]] else clust_f[i]<-max(clust_f[1:(i-1)]+1)}
table(clust_f)
boxplot(c(prob.eq))
#
baleat_f<-matrix(0,nrow=iter,ncol=length(table(clust_f))*pcov)
var_f<-matrix(0,nrow=iter,ncol=length(table(clust_f)))
for (gr in 1:length(table(clust_f))){
	ind_gr<-which(clust_f==gr)
	for (it in 1:iter){
		grupos<-S[it,ind_gr]
		var_f[it,gr]<-median(var_m[it,grupos])
		for (cov in 1:pcov){
			baleat_f[it,((gr-1)*pcov+cov)]<-mean(b_aleat[it,((grupos-1)*pcov+cov)])
}}}
#
conv<-sample(1:ncol(baleat_f),3,replace=F)
par(mfrow=c(3,1))
plot(baleat_f[,conv[1]],type='l')
plot(baleat_f[,conv[2]],type='l')
plot(baleat_f[,conv[3]],type='l')
#
par(mfrow=c(2,1))
plot(var_f[,1],type='l')
plot(var_f[,2],type='l')
#
est_Betas<-matrix(0,ncol=3,nrow=ncol(baleat_f))
for (i in 1:ncol(baleat_f)){
	est_Betas[i,2]<-mean(baleat_f[,i])
	est_Betas[i,c(1,3)]<-quantile(baleat_f[,i],c(0.025,0.975))}
#
est_var<-matrix(0,ncol=3,nrow=ncol(var_f))
for (i in 1:ncol(var_f)){
	est_var[i,2]<-median(var_f[,i])
	est_var[i,c(1,3)]<-quantile(var_f[,i],c(0.025,0.975))}
