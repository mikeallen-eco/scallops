functions {
  
  real T_dep_rec(real sbt_rec, real Topt_rec, real width_rec){
    return exp(-0.5 * ((sbt_rec - Topt_rec)/width_rec)^2); // gaussian temperature-dependent function
  }
  
  real T_dep_mort(real sbt_mort, real Topt_mort, real width_mort){
    return exp(-0.5 * ((sbt_mort - Topt_mort)/width_mort)^2); // gaussian temperature-dependent function
  }
  
  matrix age_at_length_key(real loo, real l0, real k, real cv, int n_lbins, int n_ages){
    
    vector[n_lbins] mean_age_at_length;
    vector[n_lbins] sigma_age_at_length;
    matrix[n_lbins, n_ages] prob_age_at_length;
    
    // make vector of mean ages for each length 
    for(i in 1:n_lbins){
      mean_age_at_length[i] = log((loo - i) / (loo - l0)) / -k; 
      sigma_age_at_length[i] = mean_age_at_length[i] * cv; 
    }
    
    for(i in 1:n_lbins){
      for(j in 1:n_ages){
        if(j < n_ages){
          prob_age_at_length[i,j] = normal_cdf(j+1, mean_age_at_length[i], sigma_age_at_length[i]) - normal_cdf(j, mean_age_at_length[i], sigma_age_at_length[i]);  // analog of pnorm in R
        }
        else{
          prob_age_at_length[i,j] = normal_cdf(j, mean_age_at_length[i], sigma_age_at_length[i]);
        } // close if/else
      } // close ages
    } // close lengths
    return prob_age_at_length;
  } // close function
  
} // close functions block



data {
  
  // survey data 
  
  int n_ages; // number of ages
  
  int np; // number of patches
  
  int ny_train; // years for training
  
  int ny_proj; // number of years to forecast 
  
  int n_lbins; // number of length bins (here just the range of cm values)
  
  matrix[n_ages, n_lbins] l_at_a_key;
  
  // comment if manual_selectivity == 0
  // vector[n_lbins] selectivity_at_bin;
  
  vector[n_ages] wt_at_age;
  
  real abund_p_y[np, ny_train]; // MEAN density of individuals of any age in each haul; used for rescaling the abundance to fit to our data
  
  real abund_p_y_notNA[np, ny_train]; // used for skipping NAs in abund_p_y; matrix indicating whether abund_p_y is NA (1) or not (0)

  // environmental data 
  
  real sbt_rec[np, ny_train]; // temperature data for training
  
  real sbt_mort[np, ny_train]; // temperature data for training

  real sbt_rec_proj[np, ny_proj];
  
  real sbt_mort_proj[np, ny_proj];

  real O2[np, ny_train]; // O2 data for training
  
  real O2_proj[np, ny_proj];
  
  real strat[np, ny_train]; // stratification data for training
  
  real strat_proj[np, ny_proj];
  
  // fish data
  
  real m;  // total mortality 
  
  real k;
  
  real loo;
  
  real t0;
  
  real cv;
  
  real f[n_ages, ny_train]; 
  
  real f_proj[n_ages, (ny_proj+1)];
  
  real length_50_sel_guess;
  
  vector<lower=0>[n_lbins] bin_mids;
  
  int sel_100; // age at which selectivity is 1 
  
  int age_at_maturity;
  
  int<lower = 0, upper = 1> manual_selectivity; // 0 or 1 indicating whether selectivity estimated or hard coded

  int<lower = 0, upper = 1> do_dirichlet;
  
  int<lower = 0, upper = 1> T_dep_recruitment;
  
  int<lower = 0, upper = 1> T_dep_mortality;
  
  int<lower = 0, upper = 1> eval_l_comps;
  
  int<lower = 0, upper = 1> spawner_recruit_relationship;
  
  int<lower = 0, upper = 1> run_forecast;

  int n_p_l_y[np, n_lbins, ny_train]; // SUM number of individuals in each length bin, patch, and year; used for age composition only, because the magnitude is determined by sampling effort
  
  
}

transformed data{
  
  vector[n_ages] maturity_at_age; // vector of probabilities of being mature at each age, currently binary (0/1) and taken as known
  
  for(a in 1:n_ages){
    if(a < age_at_maturity){
      maturity_at_age[a]=0;
    }else{
      maturity_at_age[a]=1;
    }
  }
  
  // matrix[n_lbins, n_ages] prob_age_at_length;
  
  // int n_p_a_y[np, n_ages, ny_train]; 
  
  // vector[n_ages] age_dist[n_lbins]; // create array of vectors 
  
  // prob_age_at_length = age_at_length_key(
    //   loo=loo,
    //   l0=l0,
    //   k=k,
    //   cv=cv,
    //   n_lbins=n_lbins,
    //   n_ages=n_ages ); 
    
    // for(p in 1:np){
      //   for(l in 1:n_lbins){
        //     for(y in 1:ny_train){
          //       age_dist[l] = prob_age_at_length[l,] * n_p_l_y[p, l, y]; // get vector of probabilities for each length, times counts
          //     }
          //   }
          // }
          
          // for(p in 1:np){
            //   for(a in 1:n_ages){
              //     for(y in 1:ny_train){
                //       n_p_a_y[p,a,y] = sum(age_dist[,a]); // not sure I'm indexing age_dist correctly, and still need to round to int somehow!
                //     }
                //   }
                // }
}

parameters{
  
  real<lower = 1e-3> sigma_r;
  
  real<lower = 1e-3> sigma_obs;
  
  real<lower=1e-3> width_rec; // sensitivity to temperature variation
  
  real Topt_rec; //  temp at which recruitment is maximized
  
  real<lower=1e-3> width_mort; // sensitivity to temperature variation
  
  real Topt_mort; //  temp at which recruitment is maximized
  
  real linbeta0rec; // intercept for logistic function vs. recruitment
  
  real O2betarec; // slope for O2 in logistic function vs. recruitment
  
  real stratbetarec; // slope for stratification in logistic function vs. recruitment

  real linbeta0mort; // intercept for logistic function vs. survival
  
  real O2betamort; // slope for O2 in logistic function vs. survival
  
  real stratbetamort; // slope for stratification in logistic function vs. survival
  
  // real<lower = -1, upper = 1> alpha; // autocorrelation term
  
  //real  log_mean_recruits; // log mean recruits per patch, changed to one value for all space/time
  vector[np] log_mean_recruits; // patch-specific mean recruits

  vector[ny_train] raw; // array of raw recruitment deviates, changed to one value per year
  
  real<upper = 0.8> p_length_50_sel; // length at 50% selectivity
  
  real<lower=0, upper=1> beta_obs; // controls how fast detection goes up with abundance
  
  // remove this if dispersal is estimated (not set to zero)
  //real<lower=0, upper=0.333> d; // dispersal fraction (0.333 = perfect admixture)
  
  real <lower = 0> theta_d;
  
  real<lower=0, upper=1> h;
  
  real log_r0;
  
}

transformed parameters{
  
  real T_adjust_rec[np, ny_train]; // tuning parameter for sbt suitability in each patch*year

  real T_adjust_mort[np, ny_train]; // tuning parameter for sbt suitability in each patch*year
  
  real length_50_sel;
  
  real sel_delta;
  
  //real mean_recruits;
  vector[np] mean_recruits;
  
  matrix<lower=0, upper=1> [np, ny_train] theta; // Bernoulli probability of encounter  
  
  real n_p_a_y_hat [np, n_ages,ny_train]; // array of numbers at patch, stage, and year 
  
  matrix[np, n_lbins] n_p_l_y_hat[ny_train]; // array number of years containing matrices with numbers at patch, length bin, and year 
  
  real dens_p_y_hat [np, ny_train]; // for tracking sum density 
  
  real dens_p_y_hat80 [np, ny_train]; // for tracking sum density >= 80 mm

  vector[ny_train-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr (it's a good or bad year everywhere)
  
  //commented while manual_selectivity==1
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  real surv[np, n_ages, ny_train];
  
  real alpha = 0;
  
  //remove this if dispersal is not set to zero
  real d = 0;
  
  real ssb0;
  
  vector[n_ages] unfished;
  
  matrix[np, ny_train] ssb;
  
  vector[n_ages] stupid_vector;
  
  real<lower=0> r0;
  
  r0 = exp(log_r0);

  ssb0 = -999;
  
  for(a in 1:n_ages){
    unfished[a] = 999;
    stupid_vector[a] = 999;
  }
  for(p in 1:np){
    for(y in 1:ny_train){
      ssb[p,y] = 999;
    }
  }
  
  if(spawner_recruit_relationship==1){
    for(a in 1:n_ages){
      if(a==1){
        unfished[a] = r0;
      }
      else{
        unfished[a] = unfished[a-1] * exp(-m);
      }
      print("unfished at age ",a," is ",unfished[a]);
      
    }
    ssb0 = sum(unfished .* maturity_at_age .* wt_at_age);
  }
  
  
  sel_delta = 2;
  
  length_50_sel = loo * p_length_50_sel; // Dan made a note to change this sometime
  
  //commented while manual_selectivity==1
  selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  // mean_selectivity_at_age = length_at_age_key * selectivity_at_bin; // calculate mean selectivity at age given variance in length at age
  
   // print("mean selectivity at age is ",mean_selectivity_at_age); // check that counts are different for every year
  
  mean_recruits = exp(log_mean_recruits);
  
   // print("mean recruits is ", mean_recruits);
   // print("sigma_obs is ", sigma_obs);
   // print("sigma_r is ", sigma_r);
  
if(T_dep_recruitment == 1){
  // calculate temperature-dependence correction factor for each patch and year depending on sbt
  for(p in 1:np){
    for(y in 1:ny_train){
      T_adjust_rec[p,y] = T_dep_rec(sbt_rec[p,y], Topt_rec, width_rec) * inv_logit(linbeta0rec + O2betarec*O2[p,y] + stratbetarec*strat[p,y]);  
      // print("T_adjust in patch", p, " and year ",y," is ",T_adjust[p,y]);

    } // close years
  } // close patches
}

if(T_dep_mortality == 1){
  // calculate temperature-dependence correction factor for each patch and year depending on sbt
  for(p in 1:np){
    for(y in 1:ny_train){
      T_adjust_mort[p,y] = T_dep_mort(sbt_mort[p,y], Topt_mort, width_mort) * inv_logit(linbeta0mort + O2betamort*O2[p,y]);  
      // print("T_adjust in patch", p, " and year ",y," is ",T_adjust[p,y]);

    } // close years
  } // close patches
}

  // print("Topt is ",Topt);
  // print("width is ", width)

  // calculate total annual mortality from instantaneous natural + fishing mortality data 
  // note that z is the proportion that survive, 1-z is the proportion that die 
  
  for(p in 1:np){
    for(a in 1:n_ages){
      for(y in 1:ny_train){
        
        if(T_dep_mortality==1){
          surv[p,a,y] = exp(-(f[a,y] + m)) * T_adjust_mort[p,y];
        }
        
        if(T_dep_mortality==0){
          surv[p,a,y] = exp(-(f[a,y] + m)) ;
        }
        
      }
    }
  }
  
  
  // fill in year 1 of n_p_a_y_hat, initialized with mean_recruits 
  for(p in 1:np){
    for(a in 1:n_ages){
      if(a==1){

        if(T_dep_recruitment==1 && spawner_recruit_relationship==0){
          //n_p_a_y_hat[p,a,1] = mean_recruits * T_adjust[p,1] * exp(raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
          n_p_a_y_hat[p,a,1] = mean_recruits[p] * T_adjust_rec[p,1] * exp(raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
        }
        if(T_dep_recruitment==0 && spawner_recruit_relationship==0){
          //n_p_a_y_hat[p,a,1] = mean_recruits * exp(raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
          n_p_a_y_hat[p,a,1] = mean_recruits[p] * exp(raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
        }
        if(T_dep_recruitment==0 && spawner_recruit_relationship==1){
          n_p_a_y_hat[p,a,1] = r0 * 0.1; // scale it down a bit -- historical fishing was still occurring
        }
        if(T_dep_recruitment==1 && spawner_recruit_relationship==1){
          n_p_a_y_hat[p,a,1] = r0 * 0.1 * T_adjust_rec[p,1];
        }
      } // close age==1 case
      else{
        n_p_a_y_hat[p,a,1] = n_p_a_y_hat[p,a-1,1] * surv[p,a-1,1]; // initialize population with mean recruitment propogated through age classes with mortality
      }
      
    } // close ages
  } // close patches
  
  
  // calculate recruitment deviates every year (not patch-specific)
  for (y in 2:ny_train){
    
    if(spawner_recruit_relationship==0){
      if (y == 2){ 
        rec_dev[y-1]  =  raw[y]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific
        // need to fix this awkward burn in
      } // close y==2 case  
      else {
        
        rec_dev[y-1] =  alpha * rec_dev[y-2] + raw[y]; // why does rec_dev[y-1] use raw[y]? 
        
      } // close ifelse
    }
    
    // describe population dynamics
    for(p in 1:np){
      
      // density-independent, temperature-dependent recruitment of age 1
      
      if(T_dep_recruitment==1 && spawner_recruit_relationship==0){
        //n_p_a_y_hat[p,1,y] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust[p,y-1];
        n_p_a_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust_rec[p,y-1];
      }
      if(T_dep_recruitment==0 && spawner_recruit_relationship==0){
        //n_p_a_y_hat[p,1,y] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) ;
        n_p_a_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev[y-1] - pow(sigma_r,2)/2) ;
      }
      
      if(T_dep_recruitment==0 && spawner_recruit_relationship==1){
        n_p_a_y_hat[p,1,y] = (0.8 * r0 * h * ssb[p, y-1]) / (0.2 * ssb0 * (1-h) + ssb0 * (h - 0.2));
      }
      if(T_dep_recruitment==1 && spawner_recruit_relationship==1){
        n_p_a_y_hat[p,1,y] = ((0.8 * r0 * h * ssb[p, y-1]) / (0.2 * ssb0 * (1-h) + ssb0 * (h - 0.2))) * T_adjust_rec[p,y-1];
      }
      // 
      // why estimate raw and sigma_r? we want to estimate process error
      // if we got rid of sigma_r, we would be saying that raw could be anything
      // that means raw could pick any value, and it would pick deviates to perfectly match abundance index
      // sigma_r scales the amount of process error that we say is possible
      // letting it go to infinity means the model will always fit the data perfectly 
      // raw is the realized process error 
      // exp(rec_dev[y-1] - pow(sigma_r,2)/2) = random variable with mean 0 and SD sigma_r
      // allows us to have a different recruitment deviation every year even though sigma_r is constant 
      
      // pop dy for non-reproductive ages 
      if(age_at_maturity > 1){ // confirm that there are non-reproductive age classes above 1
      for(a in 2:(age_at_maturity-1)){
        
        n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * surv[p,a-1,y-1]; // these just grow and die in the patch
        
      } // close ages for 2 to age at maturity
      } // close if 
      
      // pop dy for reproductive adults
      // mortality and dispersal are happening simultaneously here, between generations
      // because neither is patch-specific I don't think the order matters
      
      for(a in age_at_maturity:n_ages){
        
        // edge cases -- edges are reflecting
        if(p==1){
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * surv[p,a-1,y-1] * (1-d) + n_p_a_y_hat[p+1, a-1, y-1] * surv[p+1,a-1,y-1] * d;
        } // close patch 1 case 
        
        else if(p==np){
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * surv[p,a-1,y-1] * (1-d) + n_p_a_y_hat[p-1, a-1, y-1] * surv[p-1,a-1,y-1] * d;
        } // close highest patch
        
        else{
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * surv[p,a-1,y-1] * (1-2*d) + n_p_a_y_hat[p-1, a-1, y-1] * surv[p-1,a-1,y-1] * d + n_p_a_y_hat[p+1, a-1, y-1] * surv[p+1,a-1,y-1] * d;
          
        } // close if/else for all other patches
        
      }// close ages
    } // close patches 
    
    
  } // close year 2+ loop
  
  for(p in 1:np){
    for(y in 1:ny_train){
      
      n_p_l_y_hat[y,p,1:n_lbins] = ((l_at_a_key' * to_vector(n_p_a_y_hat[p,1:n_ages,y])) .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
      // see https://mc-stan.org/docs/2_18/reference-manual/array-data-types-section.html
      
      // n_p_l_y_hat[y,p,1:n_lbins]  =  (to_vector(n_p_l_y_hat[y,p,1:n_lbins])  .* selectivity_at_bin)';
      
      dens_p_y_hat[p,y] = sum((to_vector(n_p_l_y_hat[y,p,1:n_lbins])));
      
      dens_p_y_hat80[p,y] = sum((to_vector(n_p_l_y_hat[y,p,80:n_lbins])));

      theta[p,y] = ((1/(1+exp(-beta_obs*dens_p_y_hat[p,y]))) - 0.5)*2;
      // subtracting 0.5 and multiplying by 2 is a hacky way to get theta[0,1]
      
    }
  }
  if(spawner_recruit_relationship==1){
    for(p in 1:np){
      for(y in 1:ny_train){
        for(a in 1:n_ages){
          stupid_vector[a] = n_p_a_y_hat[p,a,y] * maturity_at_age[a] * wt_at_age[a];
        }
        ssb[p,y] = sum(stupid_vector) ;
      }
    }
  }
  
  
} // close transformed parameters block

model {
  
  real n;
  
  vector[n_lbins] prob_hat;
  
  vector[n_lbins] prob;
  
  real dml_tmp;
  
  real test;
  
  beta_obs ~ normal(0.05,0.1); // 
  
  // theta ~ uniform(0, 1); // Bernoulli probability of encounter

  // log_f ~ normal(log(m / 2),.5);
  
  if(spawner_recruit_relationship==0){
    raw ~ normal(0, sigma_r);
    sigma_r ~ normal(.7,.2); 
    
    for(i in 1:np){
    log_mean_recruits[i] ~ normal(7,3); // sd was 5
    }
  }
  
  if(spawner_recruit_relationship==1){
    h ~ normal(0.6, 0.25);
    log_r0 ~ normal(15,5);
  }
  
  Topt_rec ~ normal(9, 3);
  
  width_rec ~ normal(4, 2); 
  
  Topt_mort ~ normal(9, 3);
  
  width_mort ~ normal(4, 2); 
  
  linbeta0rec ~ normal(0, 1.4);

  O2betarec ~ normal(0, 5);
  
  stratbetarec ~ normal(0, 5);

  linbeta0mort ~ normal(0, 1.4);

  O2betamort ~ normal(0, 5);

  stratbetamort ~ normal(0, 5);
  
  // log_sigma_r ~ normal(log(.5),.1); // process error prior
  
  // alpha ~ normal(0,.25); // autocorrelation prior
  
  // uncomment this if d is estimated (not set at zero)
  //d ~ normal(0.1, 0.1); // dispersal rate as a proportion of total population size within the patch
  
  // sigma ~ normal(1,.1);  // total error prior
  // 
  // proc_ratio ~ beta(2,2);
  // 
  
  sigma_obs ~ normal(0.1,.2);  
  
  // sigma_obs ~ normal(.1, .1); 
  
  p_length_50_sel ~ normal(length_50_sel_guess/loo, .5); // sd was .2
  
  // log_scalar ~ normal(log(2),1);
  

  theta_d ~ normal(0.5,.1);
  
  
  for(y in 2:ny_train) {
    
    for(p in 1:np){
      
      if((abund_p_y[p,y]) > 0) {
        
        if(eval_l_comps==1){
          if (sum(n_p_l_y[p,1:n_lbins,y]) > 0) {
            
            if (do_dirichlet == 1){
              
              prob_hat = (to_vector(n_p_l_y_hat[y,p,1:n_lbins])  / sum(to_vector(n_p_l_y_hat[y,p,1:n_lbins])));
              
              prob = (to_vector(n_p_l_y[p,1:n_lbins,y])  / sum(to_vector(n_p_l_y[p,1:n_lbins,y])));
              
              n = sum(n_p_l_y[p,1:n_lbins,y]);
              
              dml_tmp = lgamma(n + 1) -  sum(lgamma(n * prob + 1)) + lgamma(theta_d * n) - lgamma(n + theta_d * n) + sum(lgamma(n * prob + theta_d * n * prob_hat) - lgamma(theta_d * n * prob_hat)); // see https://github.com/merrillrudd/LIME/blob/9dcfc7f7d5f56f280767c6900972de94dd1fea3b/src/LIME.cpp#L559 for log transformation of dirichlet-multinomial in Thorston et al. 2017
              
              // dml_tmp = lgamma(n + 1) - sum(lgamma(n * prob_hat + 1)) + (lgamma(theta_d * n) - lgamma(n + theta_d * n)) * prod(((lgamma(n * prob_hat + theta_d * n * prob))./(lgamma(theta_d * n * prob))));
              
              // test = prod(1:10);
              
              // print(dml_tmp);
              
              target += dml_tmp;
              
            } else {
              
              (n_p_l_y[p,1:n_lbins,y]) ~ multinomial((to_vector(n_p_l_y_hat[y,p,1:n_lbins])  / sum(to_vector(n_p_l_y_hat[y,p,1:n_lbins]))));
              
            } // close dirichlet statement
            
          } // close if any positive length comps
          
        } // close eval_length_comps
        
        log(abund_p_y[p,y]) ~ normal(log( dens_p_y_hat[p,y] + 1e-6), sigma_obs); 
        
        1 ~ bernoulli(theta[p,y]);
        
        
      } else { // only evaluate length comps if there are length comps to evaluate
      
      0 ~ bernoulli(theta[p,y]);
      
      } // close else 
    } // close patch loop
    
  } // close year loop
  
  
}


generated quantities {
  real proj_n_p_a_y_hat[np, n_ages, ny_proj+1];
  real T_adjust_rec_proj[np, ny_proj];
  real T_adjust_mort_proj[np, ny_proj];
  vector[ny_proj] rec_dev_proj;
  vector[ny_proj] raw_proj;
  real surv_proj[n_ages, (ny_proj+1)];
  matrix[np, n_lbins] proj_n_p_l_y_hat[ny_proj];
  real proj_dens_p_y_hat [np, ny_proj];
  real proj_dens_p_y_hat80 [np, ny_proj];
  // real proj_dens_p_y_hat80_lambda[np, ny_proj];
  // real dens_p_y_hat80_lambda[np, ny_train-1];

  if(run_forecast==1){
  for(p in 1:np){
    for(y in 1:ny_proj){
      T_adjust_rec_proj[p,y] = T_dep_rec(sbt_rec_proj[p,y], Topt_rec, width_rec) * inv_logit(linbeta0rec + O2betarec*O2[p,y] + stratbetarec*strat[p,y]);
    } // close years
  } // close patches
  
  for(p in 1:np){
    for(y in 1:ny_proj){
      T_adjust_mort_proj[p,y] = T_dep_mort(sbt_mort_proj[p,y], Topt_mort, width_mort) * inv_logit(linbeta0mort + O2betamort*O2[p,y]);
    } // close years
  } // close patches

  for(p in 1:np){
    for(a in 1:n_ages){
      for(y in 1:(ny_proj+1)){

        if(T_dep_mortality==0){
          surv_proj[a,y] = exp(-(f_proj[a,y] + m));
        }
        if(T_dep_mortality==1){
          	if(y==1){

	  surv_proj[a,1] = exp(-(f_proj[a,1] + m))* T_adjust_mort[p,ny_train];

	            } else {

          surv_proj[a,y] = exp(-(f_proj[a,y] + m))* T_adjust_mort_proj[p,y-1];
          
	            } // close ifelse year==1 section within T_dep_mortality==1 section

        }
      }
    }
  }

  // initialize with final year of our model
  proj_n_p_a_y_hat[,,1] = n_p_a_y_hat[,,ny_train];
  rec_dev_proj[1] = rec_dev[ny_train-1];
  raw_proj[1] = raw[ny_train];

  // project pop dy
  for(y in 2:ny_proj){
    raw_proj[y] = normal_rng(0, sigma_r);
      // print("raw_proj in year ",y," is ",raw_proj[y]);
    rec_dev_proj[y] = alpha * rec_dev_proj[y-1] + raw_proj[y];
      // print("rec_dev_proj in year ",y," is ",rec_dev_proj[y]);

  }

  for(y in 2:(ny_proj+1)){
    for(p in 1:np){

      if(T_dep_recruitment==1){
        proj_n_p_a_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2) * T_adjust_rec_proj[p,y-1];
      }
      if(T_dep_recruitment==0){
        proj_n_p_a_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2);
      }

      if(age_at_maturity > 1){
        for(a in 2:(age_at_maturity-1)){
          proj_n_p_a_y_hat[p,a,y] = proj_n_p_a_y_hat[p, a-1, y-1] * surv_proj[a-1,y-1];
        } // close ages for 2 to age at maturity
      } // close if

      for(a in age_at_maturity:n_ages){
        if(p==1){
          proj_n_p_a_y_hat[p,a,y] = proj_n_p_a_y_hat[p, a-1, y-1] * surv_proj[a-1,y-1] * (1-d) + proj_n_p_a_y_hat[p+1, a-1, y-1] * surv_proj[a-1,y-1] * d;
        } // close patch 1 case

        else if(p==np){
          proj_n_p_a_y_hat[p,a,y] = proj_n_p_a_y_hat[p, a-1, y-1] * surv_proj[a-1,y-1] * (1-d) + proj_n_p_a_y_hat[p-1, a-1, y-1] * surv_proj[a-1,y-1] * d;
        } // close highest patch

        else{
          proj_n_p_a_y_hat[p,a,y] = proj_n_p_a_y_hat[p, a-1, y-1] * surv_proj[a-1,y-1] * (1-2*d) + proj_n_p_a_y_hat[p-1, a-1, y-1] * surv_proj[a-1,y-1] * d + proj_n_p_a_y_hat[p+1, a-1, y-1] * surv_proj[a-1,y-1] * d;

        } // close if/else for all other patches

      }// close ages
    } // close patches


  } // close year 2+ loop

  for(p in 1:np){
    for(y in 1:(ny_proj)){

      proj_n_p_l_y_hat[y,p,1:n_lbins] = ((l_at_a_key' * to_vector(proj_n_p_a_y_hat[p,1:n_ages,y])) .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
      proj_dens_p_y_hat[p,y] = sum((to_vector(proj_n_p_l_y_hat[y,p,1:n_lbins])));
      proj_dens_p_y_hat80[p,y] = sum((to_vector(proj_n_p_l_y_hat[y,p,80:n_lbins])));
    
    }
  }
  
  // // to estimate projected annual population change
  // for(p in 1:np){
  //   
  // proj_dens_p_y_hat80_lambda[p,1] = proj_dens_p_y_hat80[p,1] / dens_p_y_hat[p,ny_train];
  // 
  // for(y in 2:(ny_proj)){
  // 
  // proj_dens_p_y_hat80_lambda[p,y] = proj_dens_p_y_hat80[p,y] / proj_dens_p_y_hat80[p,y-1];
  // 
  //   }
  // }
  // 
  // // to estimate observed annual population change
  // for(p in 1:np){
  //   for(y in 1:(ny_train-1)){
  // 
  // dens_p_y_hat80_lambda[p,y] = dens_p_y_hat80[p,y+1] / dens_p_y_hat80[p,y];
  // 
  //   }
  // }
    
}
}

