#Install Rstan
setwd("../")
library(rstan)
library(MCMCvis)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

cam_allpH_data <-"
data {
	int N;    
	int K;
	vector[N] logNt;     
	real t[N];
	real<lower=3.0, upper=3.8> pH[N];    
	int<lower=1, upper=K> strain_index[N];    
	real Nu;
}	

parameters {
	vector<lower=4, upper=7.5>[K] logN0;
	real<lower=4, upper=7.5> mean_logN0;  
	real<lower=0> sigma_logN0;
	real<lower=0, upper=1> sigma_logNt;
	vector[4] abef[K];
	vector[4] abef0;
	cholesky_factor_corr[4] corr_chol;
	vector<lower=0>[4] sigma_vec;
}
transformed parameters{
	vector[K] a;
	vector[K] b;
	vector[K] e;
	vector[K] f;
	cholesky_factor_cov[4] cov_chol;
	vector[N] mean_survivor;
	vector[N] d;
	vector[N] p;
	for(k in 1:K){
		a[k] = abef[k,1];
		b[k] = abef[k,2];
		e[k] = abef[k,3];
		f[k] = abef[k,4];
	}
	cov_chol = diag_pre_multiply(sigma_vec,corr_chol);
	for(n in 1:N)
	{
		d[n] = exp(a[strain_index[n]]*pH[n]+b[strain_index[n]]);
		p[n] = exp(e[strain_index[n]]*pH[n]+f[strain_index[n]]);
		mean_survivor[n] = logN0[strain_index[n]] -(t[n]/d[n])^p[n];
	}
}
model {
	logN0 ~ normal(mean_logN0, sigma_logN0);
	corr_chol ~ lkj_corr_cholesky(Nu);
	abef ~ multi_normal_cholesky(abef0, cov_chol);
	logNt ~ normal(mean_survivor, sigma_logNt);
}
generated quantities {
	matrix[4,4] corr;
	matrix[4,4] cov;
	vector[N] log_lik;
	corr = multiply_lower_tri_self_transpose(corr_chol);
	cov = multiply_lower_tri_self_transpose(cov_chol);
	for(n in 1:N){
		log_lik[n] = normal_lpdf(mean_survivor[n]|logN0[strain_index[n]] -(t[n]/d[n])^p[n],sigma_logNt);
	}
}
"

#Load data
d <- read.csv("data/All_dataset0.csv")
d <- subset(d, d$pH != 3.3)
d <- subset(d, d$pH != 3.5)

#Convert data to list
listd <- list(t = as.numeric(d$time), 
              logNt = as.numeric(d$logN),
              strain_index = as.numeric(d$strain),
              pH = as.numeric(d$pH), 
              N = length(d$time), 
              K = 11, 
              Nu = 1
              )

#MCMC_fitting
cam_allpH_stan <- stan( model_code = cam_allpH_data,
                        data = listd,
                        iter = 5000,
                        warmup = 2500,
                        chain = 4,                  
                        thin = 1,
                        seed=123,
                        control = list(adapt_delta = 0.8, max_treedepth = 15)
                        )

#MCMC_output
ms <- rstan::extract(cam_allpH_stan)
saveRDS(ms, file = 'result_Model_para.obj')

write.csv(MCMCsummary(cam_allpH_stan),"MCMC_summary_para.csv")