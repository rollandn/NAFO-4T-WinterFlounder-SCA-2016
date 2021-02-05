solve_Baranov <- function(Bprime, S, C, M, nIter){
   # Solve Baranov catch equation:
         
   F = C / sum(Bprime); # Initial estimate of 'F'.
   ZNew = rep(NA, n_age);
   for (i in 1:nIter){
      for (a in 1:n_age) ZNew[a] = M[a] + S[a] * F;  # Vector.
     
      f = 0;
      for (a in 1:n_age) f = f + (Bprime[a] * F) * (1 - exp(-ZNew[a])) / ZNew[a];
      f = f - C;    
     
      J = 0;
      for (a in 1:n_age) J = J + (Bprime[a] / ZNew[a]) * ((S[a] * ZNew[a] * F) + (1 - exp(-ZNew[a])) - ((S[a] * F / ZNew[a]) * (1-exp(-ZNew[a]))));    
      
      F = F - f / J;
   }
   
   return(F)
}
