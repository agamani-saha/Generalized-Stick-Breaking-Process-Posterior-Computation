###########################################################################################
library(MCMCprecision)
library(markovchain)
library(ggplot2)
set.seed(1)
d = 15
c = 9
########################### True TPM ######################################################
p = rep(NA, d-1)
trans.matrix.true = matrix(rep(NA, d*d), nrow = d)
for(i in 0:(d-1)){
  for(j in 0:(d-1)){
    p = 1/(log(i+c))
    trans.matrix.true[i+1,j+1] = ((1-p)^j)*p
  }
}
trans.matrix.true
#############################################################################################################
###################################### GET DATASET ######################################################################
# Install and load GSODR (only once)
library(bit64)
library(dplyr)
library(lubridate)
library(MCMCprecision)
library(markovchain)
library(GSODR)

heathrow <- read.csv("HeathrowRainMonthly1995To31Aug2025.csv")
heathrow_rain <- heathrow[, c("YearMonth","PRCP","states")]

###################################### FUNCTION TO COMPUTE TPM ##############################################
compute_tpm <- function(states, prob) {
  # uniq_states <- 0:max(monthly_rain$states)
  uniq_states <- 0:max(heathrow_rain$states)
  d <- length(uniq_states)
  
  # Initialize transition matrix
  tpm <- matrix(0, nrow = d, ncol = d, 
                dimnames = list(uniq_states, uniq_states))
  
  # Count transitions
  for (i in 1:(length(states) - 1)) {
    tpm[as.character(states[i]), as.character(states[i + 1])] <-
      tpm[as.character(states[i]), as.character(states[i + 1])] + 1
  }
  
  # Convert counts to probabilities
  if(prob == TRUE){
    tpm <- tpm / rowSums(tpm)
  } 
  tpm[is.na(tpm)] <- 0
  return(tpm)
}
###################################### COMPUTE TPM ##############################################
TPM_full <- compute_tpm(heathrow_rain$states, TRUE)
###################################### OUTPUT TPM ##############################################
cat("\n--- Full Dataset TPM ---\n")
print(round(TPM_full,3))
################################################################################################
p_hat <- rep(NA,d)
for(i in 1:(d+1)){
  df <- data.frame(value = 0:(ncol(TPM_full)-1), prob = TPM_full[i,])
  # empirical mean
  mu_hat <- sum(df$value * df$prob)
  # fitted p (support starts at 0)
  p_hat[i] <- 1 / (1 + mu_hat)
}
matplot(1:(d+1), cbind(p_hat, p), type = "l")
i = 1:(d+1)
x = 1/log(i+c)
fit <- lm(p_hat ~ x-1)
summary(fit)
################################## Improved scientific-style plot ###############
library(ggplot2)
library(patchwork)
library(latex2exp)
library(reshape2)

## ---------- Data for Plot 1: Fitted vs Estimated p ----------
df1 <- data.frame(
  State = 1:length(fitted(fit)),
  Fitted = fitted(fit),
  Estimated = p_hat
)

df1_long <- melt(df1, id.vars = "State",
                 variable.name = "Type", value.name = "Probability")

p1 <- ggplot(df1_long, aes(x = State, y = Probability, 
                           color = Type, linetype = Type)) +
  geom_line(linewidth = 0.7) +   # thicker lines for visibility
  scale_color_manual(values = c("Fitted" = "brown", "Estimated" = "blue")) +
  scale_linetype_manual(values = c("Fitted" = "solid", "Estimated" = "dashed")) +
  labs(x = "States", y = TeX("Probability ( $p$, $\\hat{p}$ )")) +
  theme_minimal(base_size = 16) +   # <-- controls overall text size
  theme(
    legend.position = c(0.8, 0.85),
    legend.background = element_blank(),
    legend.title = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

## ---------- Data for Plot 2: Geometric Fit ----------
df2 <- data.frame(value = 0:(ncol(TPM_full)-1), prob = TPM_full[3, ])
mu_hat <- sum(df2$value * df2$prob)
p_hat <- 1 / (1 + mu_hat)
df2$fitted_prob <- dgeom(df2$value, prob = p_hat)

df2_long <- melt(df2[, c("value","prob","fitted_prob")], id.vars = "value",
                 variable.name = "Type", value.name = "Probability")

df2_long$Type <- factor(df2_long$Type,
                        levels = c("prob","fitted_prob"),
                        labels = c("Actual","Fitted"))

p2 <- ggplot(df2_long, aes(x = value, y = Probability, 
                           color = Type, linetype = Type)) +
  geom_line(linewidth = 0.7) +   # thicker lines
  scale_color_manual(values = c("Actual" = "brown", "Fitted" = "blue")) +
  scale_linetype_manual(values = c("Actual" = "solid", "Fitted" = "dashed")) +
  labs(x = "States", y = "Probability") +
  theme_minimal(base_size = 16) +   # <-- overall text size
  theme(
    legend.position = c(0.8, 0.85),
    legend.background = element_blank(),
    legend.title = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

## ---------- Combine Side by Side ----------
final_plot <- p1 | p2
final_plot

