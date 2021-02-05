// Yellowtail stage-structured population model:
#include <TMB.hpp>
   
template<class Type>
Type objective_function<Type>::operator()(){
   // ============================== DATA SECTION ==============================
   // Variable declarations:   
   DATA_VECTOR(C);                     // Landings in tonnes (n_year).       
   DATA_VECTOR(B_rv);                  // Biomass index from survey (in MT). 
   DATA_MATRIX(P_rv);                  // Survey proportions-at-age (n_year x n_age) (by number).
   DATA_MATRIX(P_fishery);             // Fishery proportions-at-age (n_year x n_age) (by number) 'LgrpObsProp'.
   DATA_MATRIX(weight_at_age_rv);      // Observed mean weights at age by year (n_year x n_age).  
   DATA_MATRIX(weight_at_age_fishery); // Observed mean weights at age by year (n_year x n_age).  
   DATA_MATRIX(maturity_at_age);       // observed proportions mature by age and year (n_year x n_age).   
   DATA_VECTOR(ages);                  // Age vector (e.g. 1...8).
   DATA_FACTOR(S_rv_block);            // Survey selectivity block mapping index.
   DATA_FACTOR(S_fishery_block);       // Fishery selectivity block mapping index.
   DATA_FACTOR(M_block_year);          // Mortality block mapping index by year.
   DATA_FACTOR(M_block_age);           // Mortality block mapping index by age.
   // DATA_INTEGER(splitAge);          // Index of split age in the 'ages' vector.
   DATA_SCALAR(frac_rv);               // Fraction of year where the RV survey occurs.
   DATA_INTEGER(baranovIter);          // Number of Baranov iterations. 
   // DATA_VECTOR(loglike_weight);     // Log-likelihood component weightings (RVjuv, RVadt, CatchLgrps).
   // DATA_SCALAR(init_tauLgrp);       // Initial values for Lgrp prop observation std errors.

   // ========================== PARAMETER_SECTION =============================
   // Recruitment parameters:
   PARAMETER(log_AvgR);                // Log-scale average recruitment:
   Type avgR = exp(log_AvgR);          // Average recruitment.
   PARAMETER(sigma_R);                 // 'R_priorSD' recruitment (= 0.5).
   PARAMETER_VECTOR(init_recDevs);     // ages [-5,5], (sage,lage) Recruitment deviations: initial abundance 1985    
   PARAMETER_VECTOR(recDevs);          // Recruitment deviations: 1977-2014. (2 x nT)     
   PARAMETER(logit_gamma_R);           // Logit-scale recruitment correlation parameter.
   Type gamma_R = exp(logit_gamma_R) / (Type(1) + exp(logit_gamma_R)); // Recruitment correlation parameter.

   // Fishery mortality parameters:
   PARAMETER(log_Finit);               // Log-scale initial fishing mortality prior.   
   Type Finit = exp(log_Finit);        // Initial fishing mortality prior.
   
   // Survey catchability parmeters:
   PARAMETER(log_q_rv);                // Log-scale survey catchability coefficient 'rvq'. (-5, 1.5, ph_q);
   Type q_rv = exp(log_q_rv);          // Survey catchability coefficient 'rvq'.
   PARAMETER(log_q_prior_rv);          // Survey catchability prior.
   PARAMETER(sigma_q_rv);              // Survey catchability error 0.15.
   
   // Natural morality parameters:
   PARAMETER_VECTOR(Minit_prior);      // natural mortality Initial values (1x2) for ages 1-3, 5-8. 
   PARAMETER_VECTOR(sigma_M);          // 'M_priorSD' Prior over standard errors for M (1x2);
   PARAMETER_MATRIX(log_M);            // Log-scale natural moraltity parameters.
 
   // Selectivity curve parameters:
   PARAMETER_VECTOR(log_S50);          // Log-scale logistic selectivity parameter (age at 50% level). 
   PARAMETER_VECTOR(log_S95_delta);    // Log-scale logistic selectivity delta parameter (defines age at 5% and 95% level). 

   // Survey index parameters:
   PARAMETER(log_sigma_rv);           
   Type sigma_rv = exp(log_sigma_rv);  // Survey error parameter.

   // Define data dimensions:
   int n_year = C.size();              // Number of years in time series.
   int n_age = ages.size();            // Number of age categories.
 
   // Mortality mapping code block:   
   matrix<Type> M(n_year, n_age);      // Natural mortality matrix.
   matrix<Type> Z(n_year, n_age);      // Total mortality matrix.
   for (int t = 0; t < n_year; t++){
      for (int a = 0; a < n_age; a++){
         M(t,a) = exp(log_M(M_block_year[t]-1, M_block_age[a]-1)); // Map natural mortality parameters.
         Z(t,a) = M(t,a); // Initial total mortality set to natural mortality. 
      }
   } 
   
   // Selectivity mapping code block: 
   matrix<Type> S_rv(n_year, n_age); 
   matrix<Type> S_fishery(n_year, n_age);   
   vector<Type> S50_rv(n_year);
   vector<Type> S50_fishery(n_year); 
   vector<Type> S95_rv(n_year);
   vector<Type> S95_fishery(n_year); 
 
   // Map selectivity curve parameters to proper survey and fishery block:
   for (int i = 0; i < n_year; i++){
      S50_rv[i] = exp(log_S50[S_rv_block[i]-1]);
      S95_rv[i] = S50_rv[i] + exp(log_S95_delta[S_rv_block[i]-1]);
      S50_fishery[i] = exp(log_S50[S_fishery_block[i]-1]);
      S95_fishery[i] = S50_fishery[i] + exp(log_S95_delta[S_fishery_block[i]-1]);        
   }
      
   // Define full selectivity matrices (year x age):
   for (int i = 0; i < n_year; i++){
      Type tmp_rv = log(Type(19)) / (S95_rv[i] - S50_rv[i]);
      Type tmp_fishery = log(Type(19)) / (S95_fishery[i] - S50_fishery[i]); 
      for (int j = 0; j < n_age; j++){
         S_rv(i,j) = Type(1) / (Type(1) + exp(-tmp_rv * (ages[j] - S50_rv[i])));   
         S_fishery(i,j) = Type(1) / (Type(1) + exp(-tmp_fishery * (ages[j] - S50_fishery[i])));
      }
   }

   // Calculate population abudance and biomass by length group for first year:
   matrix<Type> N(n_year, n_age);              // Population numbers (in thousands).
   matrix<Type> N_rv(n_year, n_age);           // Survey trawlable abundance (in thousands).
   matrix<Type> B_rv_pred(n_year, n_age);      // Survey trawlable biomass (MT).
   vector<Type> B_rv_pred_year(n_year);        // Total annual survey trawlable biomass (MT).
   matrix<Type> N_fishery(n_year, n_age);      // Fishable abundance (in thousands).
   matrix<Type> SSB(n_year, n_age);            // Spawning stock biomass.
   vector<Type> SSB_year(n_year);              // Total annual spawning stock biomass.
   vector<Type> Bprime(n_age);                 // Fishable biomass for Baranov iterations.
   matrix<Type> P_rv_pred(n_year, n_age);      // Predicted proportion-at-age in survey.
   matrix<Type> P_fishery_pred(n_year, n_age); // Predicted proportion-at-age in fishery ('uCgat').
   vector<Type> F(n_year);                     // Fishing mortality.
      
   //================= Population dynamics code block ==========================
   // Calculate poplation numbers for first year:
   N(0,0) = avgR * exp(init_recDevs[0]); // Recruitment deviation for first year.
   for (int a = 1; a < n_age; a++){
      N(0,a) = log(avgR) + init_recDevs[a]; // Independent recruitment deviations
      for (int j = 0; j < a; j++) N(0,a) -= S_fishery(0,j) * Finit + M(0,j);  
      N(0,a) = exp(N(0,a));
   }
   N(0,n_age-1) = N(0,n_age-1) / (Type(1.0) - exp(-(S_fishery(0,n_age-1) * Finit + M(0,n_age-1))));  // Last age.
   
   // Calculate exploitable biomass 'Bprime' for first year:
   for (int a = 0; a < n_age; a++) Bprime[a] = S_fishery(0,a) * N(0,a) * weight_at_age_fishery(0,a);  
   
   // Solve Baranov catch equation for first year: 
   F[0] = C[0] / sum(Bprime); // Initial estimate of 'F'.
   vector<Type> ZNew(n_age);
   for (int i = 0; i < baranovIter; i++){
      Type f = 0; Type J = 0;
      for (int a = 0; a < n_age; a++) ZNew[a] = M(0,a) + S_fishery(0,a) * F[0];  
      for (int a = 0; a < n_age; a++) f += (Bprime[a] * F[0]) * (1 - exp(-ZNew[a])) / ZNew[a];
      f -= C[0];    
      for (int a = 0; a < n_age; a++)
         J += (Bprime[a] / ZNew[a]) * ((S_fishery(0,a) * ZNew[a] * F[0]) + (Type(1) - exp(-ZNew[a])) - ((S_fishery(0,a) * F[0] / ZNew[a]) * (Type(1)-exp(-ZNew[a]))));   
      F[0] -= f / J;
   }
   for (int a = 0; a < n_age; a++) Z(0,a) =  M(0,a) + S_fishery(0,a) * F[0]; // Update total mortality.
      
   // Loop over remaining years:
   for (int t = 1; t < n_year; t++){
      // Calculate population numbers for subsequent years:
      N(t,0) = exp(gamma_R * log(N(t-1,0)) + (Type(1)-gamma_R) * log_AvgR + recDevs(t-1));          // Recruitment for first age.
      for (int a = 1; a < (n_age-1); a++ ) N(t,a) = N(t-1,a-1) * exp(-Z(t-1,a-1));                  // Abundance intermediate ages.
      N(t,n_age-1) = N(t-1,n_age-2) * exp(-Z(t-1,n_age-2)) + N(t-1,n_age-1) * exp(-Z(t-1,n_age-1)); // Abundance last age. 
      
      // Calculate exploitable biomass 'Bprime':
      for (int a = 0; a < n_age; a++) Bprime[a] = S_fishery(t,a) * N(t,a) * weight_at_age_fishery(t,a);  
      
      // Solve Baranov catch equation for each year: 
      F[t] = C[t] / sum(Bprime); // Initial estimate of 'F'.
      vector<Type> ZNew(n_age);
      for (int i = 0; i < baranovIter; i++){
         Type f = 0; Type J = 0;
         for (int a = 0; a < n_age; a++) ZNew[a] = M(t,a) + S_fishery(t,a) * F[t];  // Vector.
         for (int a = 0; a < n_age; a++) f += (Bprime[a] * F[t]) * (1 - exp(-ZNew[a])) / ZNew[a];
         f -= C[t];    
         for (int a = 0; a < n_age; a++)
            J += (Bprime[a] / ZNew[a]) * ((S_fishery(t,a) * ZNew[a] * F[t]) + (Type(1) - exp(-ZNew[a])) - ((S_fishery(t,a) * F[t] / ZNew[a]) * (Type(1)-exp(-ZNew[a]))));   
         F[t] -= f / J;
      }
      for (int a = 0; a < n_age; a++) Z(t,a) =  M(t,a) + S_fishery(t,a) * F[t]; // Update total mortality.
   }  

   // Calculate summary abundance and biomass values:
   B_rv_pred.fill(0);
   B_rv_pred_year.fill(0);
   SSB_year.fill(0);
   for (int t = 0; t < n_year; t++){
      for (int a = 0; a < n_age; a++){
         // Calculate Spawning Stock Biomass (SSB) values:
         SSB(t,a) = maturity_at_age(t,a) * N(t,a) * weight_at_age_rv(t,a);
         SSB_year[t] += SSB(t,a); 
         
         // Calculate predicted RV survey abundance:
         N_rv(t,a) = q_rv * S_rv(t,a) * N(t,a) * exp(-frac_rv * Z(t,a)); 
         
         // Survey trawlable biomass (MT).
         B_rv_pred(t,a) = N_rv(t,a) * weight_at_age_rv(t,a);           
 
         // Calculate total annual predicted RV survey total biomass (MT):
         B_rv_pred_year[t] += B_rv_pred(t,a); 
         
         // Calculate predicted fishable abundance:
         N_fishery(t,a) = S_fishery(t,a) * N(t,a);
      }     
   }
   
   // Calculate predicted survey catch proportions and RV residuals:
   vector<Type> N_rv_year(n_year);
   N_rv_year.fill(0);
   for (int t = 0; t < n_year; t++){
      for (int a = 0; a < n_age; a++) N_rv_year[t] += N_rv(t,a);
      for (int a = 0; a < n_age; a++) P_rv_pred(t,a) = N_rv(t,a) / N_rv_year[t];                                              
   }
   
   // Calculate predicted fishery catch proportions:
   vector<Type> N_fishery_year(n_year);
   N_fishery_year.fill(0);
   for (int t = 0; t < n_year; t++){
      for (int a = 0; a < n_age; a++) N_fishery_year[t] += N_fishery(t,a);
      for (int a = 0; a < n_age; a++) P_fishery_pred(t,a) = N_fishery(t,a) / N_fishery_year[t];
   }
                  
   //====================== Calculate Index Likelihood =========================
   vector<Type> rv_residuals_year(n_year);   
   Type loglike_index = 0;
   for (int t = 0; t < n_year; t++){
      rv_residuals_year[t] = log(B_rv[t]) - log(B_rv_pred_year[t]);
      loglike_index += log(sigma_rv) + (pow(rv_residuals_year[t], 2) / (Type(2.0) * pow(sigma_rv,2))); 
   } 
   
   //============== Calculate RV Index Age Proportion Likelihood ================
   Type loglike_P_rv = 0;
   Type loglike_P_rv_zero = 0;
   vector<Type> res(n_age);
   Type etaSumSq = 0;
   for (int t = 0; t < n_year; t++){           
      int n = 0;
      res.fill(0);
      for (int a = 0; a < n_age; a++){
         if (P_rv(t,a) > 0.00001){
            res[a] = log(P_rv(t,a)) - log(P_rv_pred(t,a));
            n++; 
         }else{
            // Large penalty for zero proportions:
            loglike_P_rv_zero += -Type(50) * log(Type(1) - P_rv_pred(t,a)); 
         }
      } 
      Type mu_res = sum(res) / n;
      for (int a = 0; a < n_age; a++){
         if (P_rv(t,a) > 0.00001) etaSumSq += pow(res[a] - mu_res, 2);
      }
   }  
   loglike_P_rv += n_year * (log(etaSumSq) - log(n_year));
   
   //============== Calculate Fishery Age Proportion Likelihood ================
   Type loglike_P_fishery = 0;
   Type loglike_P_fishery_zero = 0;
   etaSumSq = 0;
   for (int t = 0; t < n_year; t++){           
      int n = 0;
      res.fill(0);
      for (int a = 0; a < n_age; a++){
         if (P_fishery(t,a) > 0.00001){
            res[a] = log(P_fishery(t,a)) - log(P_fishery_pred(t,a));
            n++; 
         }else{
            // Large penalty for zero proportions:
            loglike_P_fishery_zero += -Type(50) * log(Type(1) - P_fishery_pred(t,a)); 
         }
      } 
      Type mu_res = sum(res) / n;
      for (int a = 0; a < n_age; a++){
         if (P_fishery(t,a) > 0.00001) etaSumSq += pow(res[a] - mu_res, 2);
      }
   }  
   loglike_P_fishery += n_year * (log(etaSumSq) - log(n_year));
   
   //================= Calculate Recruitment Prior Likelihood ==================
   // Initialize recruitment likelihood:
   Type loglike_R = 0;
   for (int j = 0; j < n_age; j++) loglike_R += log(sigma_R) + 0.5 * (init_recDevs[j] * init_recDevs[j]) / (sigma_R * sigma_R); // first year.
   for (int i = 0; i < (n_year-1); i++) loglike_R += log(sigma_R) + 0.5 * (recDevs[i] * recDevs[i]) / (sigma_R * sigma_R);      // first age.

   //================== Calculate Mortality Prior Likelihood ===================
   // Normal prior for initial M in first block:
   Type loglike_M = 0;
   for (int j = 0; j < log_M.cols(); j++) loglike_M += log(sigma_M[j]) + 0.5 * pow(exp(log_M(0,j)) - Minit_prior[j],2) / (sigma_M[j] * sigma_M[j]);
  
   //================ Calculate Catchability Prior Likelihood ==================
   // Survey catchability 'q' prior:
   Type loglike_q = 0;
   loglike_q  += 0.5 * pow(log_q_rv - log_q_prior_rv, 2) / (sigma_q_rv * sigma_q_rv);
   
   //====================== Calculate Exploitation Rate ========================
   // Revised code:
   //matrix<Type> subHt(n_year, n_age);
   //vector<Type> p(2);
   //for (int t = 0; t < n_year; t++){
   //   for (int g = 0; g < 2; g++) p[g] = LgrpObsProp(t,g) * Cw(t,g);
   //   p = p / sum(p);
   //   for (int g = 0; g < 2; g++) subHt(t,g) = (p[g] * C[t]) / bLgrp(t,g); // Harvest rate matrix.
   //}  

   // Output model variables:
   REPORT(rv_residuals_year);
   REPORT(N);
   REPORT(N_rv);  
   REPORT(N_fishery);
   REPORT(B_rv_pred);
   REPORT(B_rv_pred_year);
   REPORT(P_rv_pred);
   REPORT(P_fishery_pred);
   REPORT(S_rv);
   REPORT(S_fishery);
   REPORT(M);
   REPORT(F);
   REPORT(SSB);
   REPORT(SSB_year);
   REPORT(P_fishery_pred);
   REPORT(loglike_index);
   REPORT(loglike_P_rv);
   REPORT(loglike_P_rv_zero);
   REPORT(loglike_P_fishery);
   REPORT(loglike_P_fishery_zero);
   REPORT(loglike_R);
   REPORT(loglike_M);
   REPORT(loglike_q);
   
   // Loglikelihood contributions: 
   Type loglike = 0;
   loglike += loglike_index; 
   loglike += loglike_P_rv;
   loglike += loglike_P_rv_zero;
   loglike += loglike_P_fishery;
   loglike += loglike_P_fishery_zero;
   loglike += loglike_R; 
   loglike += loglike_M;
   loglike += loglike_q;
   
   return loglike;
}

