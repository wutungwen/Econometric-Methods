#Implement the method introduced in Imbens, G. W. and D. B. Rubin. "Estimating Outcome Distributions for Compliers in Instrumental Variables Models." The Review of Economic Studies 64, no. 4 (1997): 555-574. and visualization of Abadie, Alberto. "Bootstrap Tests for Distributional Treatment Effects in Instrumental Variable Models." Journal of the American Statistical Association 97, no. 457 (2002): 284-292.

#the column names should be "outcome", "treatment"(binary), and "instrument"(binary)!

outcome_distributions <- function(data, binwidth){
  Y1_bar <- mean(subset(data, instrument == 1)$outcome)
  Y0_bar <- mean(subset(data, instrument == 0)$outcome)
  D1_bar <- mean(subset(data, instrument == 1)$treatment)
  D0_bar <- mean(subset(data, instrument == 0)$treatment)
  
  f00_df <- subset(data, instrument == 0 & treatment == 0)
  f01_df <- subset(data, instrument == 0 & treatment == 1)
  f10_df <- subset(data, instrument == 1 & treatment == 0)
  f11_df <- subset(data, instrument == 1 & treatment == 1)
  
  phi_n <- 1 - mean(subset(data, instrument == 1)$treatment)
  phi_a <- mean(subset(data, instrument == 0)$treatment)
  phi_c <- 1 - phi_n - phi_a
  
  bin <- seq(min(data$outcome), max(data$outcome), by = binwidth)
  F00 <- c()
  F01 <- c()
  F10 <- c()
  F11 <- c()
  
  for (i in 1:length(bin)){
    F00[i] <- nrow(subset(f00_df, outcome < bin[i]))
    F01[i] <- nrow(subset(f01_df, outcome < bin[i]))
    F10[i] <- nrow(subset(f10_df, outcome < bin[i]))
    F11[i] <- nrow(subset(f11_df, outcome < bin[i]))
  }
  
  g_c0 <- c()
  g_c1 <- c()
  
  for (i in 1:length(bin)){
    g_c0[i] <- ((phi_n + phi_c)/(phi_c))*F00[i] - (phi_n/phi_c)*F10[i]
    g_c1[i] <- ((phi_a + phi_c)/(phi_c))*F11[i] - (phi_a/phi_c)*F01[i]
  }
  
  g_c0_scaled <- c()
  g_c1_scaled <- c()
  
  for (i in 1:length(g_c0)){
    g_c0_scaled[i] <- g_c0[i]/max(g_c0)
    g_c1_scaled[i] <- g_c1[i]/max(g_c1)
  }
  return (list(Y1_bar = Y1_bar, Y0_bar = Y0_bar, D1_bar = D1_bar, D0_bar = D0_bar,
               phi_n = phi_n, phi_a = phi_a, phi_c = phi_c,
               F00 = F00, F10 = F10, F01 = F01, F11 = F11,
               g_c0 = g_c0_scaled, g_c1 = g_c1_scaled, bin = bin))
}
