library(quantmod)
c = 10
###################################### GET DATASET ##############################################
getSymbols("USO", src = "yahoo", from = "2015-01-01", to = "2020-01-01")
volumes <- Vo(USO)

# Convert to data frame and divide by 1e6
vol_data <- data.frame(Date = index(volumes), Volume = as.numeric(volumes) / 1e6)

# Fix date display
start_date <- format(as.Date(index(volumes)[1]), "%Y-%m-%d")
end_date <- format(as.Date(index(volumes)[nrow(volumes)]), "%Y-%m-%d")
cat("USO data available from", start_date, "to", end_date, "\n")

###################################### CREATE STATES (Floor Integers) ##############################################
vol_data$State <- floor(vol_data$Volume)

###################################### FUNCTION TO COMPUTE TPM ##############################################
compute_tpm <- function(states, prob) {
  d <- max(states) + 1
  uniq_states <- 0:(d-1)
  
  # Initialize transition matrix
  tpm <- matrix(0, nrow = d, ncol = d, dimnames = list(uniq_states, uniq_states))
  
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
TPM_full <- compute_tpm(vol_data$State, TRUE)
###################################### OUTPUT TPM ##############################################
cat("\n--- Full Dataset TPM ---\n")
print(TPM_full)
########################################### TRUE TPM #############################################
d = nrow(TPM_full)
trans.matrix.true = matrix(rep(NA, d*d), nrow = d)
for(i in 0:(d-1)){
  for(j in 0:(d-1)){
    lambda = (log(i+c))
    trans.matrix.true[i+1,j+1] = exp(-lambda)*(lambda^j)/factorial(j)
  }
}
trans.matrix.true
###################################### PLOT LAMBDA GRAPH ##############################################
i = 0:(d-1)
lambda = log(i+c)
plot(lambda)
################################# 
lambda_hat <- rep(NA,d)
for(i in 1:d){
  n = 0:(d-1)
  lambda_hat[i] <- sum(n*TPM_full[i,])
}
lambda_hat
plot(lambda_hat)
matplot(1:d, cbind(lambda_hat, lambda[1:14]), type = "l")
i = 1:14
x = log(i+c)
fit <- lm(lambda_hat ~ x)
summary(fit)
############################# GRAPH ################################################
library(ggplot2)
library(patchwork)
library(latex2exp)
library(reshape2)

## ---------- Data for Plot 1: Fitted λ vs Estimated λ ----------
df1 <- data.frame(
  State = 1:length(fitted(fit)),
  Fitted = fitted(fit),
  Estimated = lambda_hat[1:14]
)

df1_long <- melt(df1, id.vars = "State",
                 variable.name = "Type", value.name = "Rate")

p1 <- ggplot(df1_long, aes(x = State, y = Rate, 
                           color = Type, linetype = Type)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("Fitted" = "brown", "Estimated" = "blue")) +
  scale_linetype_manual(values = c("Fitted" = "solid", "Estimated" = "dashed")) +
  labs(x = "States", y = TeX("Rate ($\\lambda$, $\\hat{\\lambda}$)")) +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(5, 5, 5, 5))

## ---------- Data for Plot 2: Poisson Fit ----------
df2 <- data.frame(value = 0:(ncol(TPM_full)-1), prob = TPM_full[2,])
lambda_hat_val <- sum(df2$value * df2$prob)
df2$fitted_prob <- dpois(df2$value, lambda = lambda_hat_val)

df2_long <- melt(df2[, c("value","prob","fitted_prob")], id.vars = "value",
                 variable.name = "Type", value.name = "Probability")

df2_long$Type <- factor(df2_long$Type,
                        levels = c("prob","fitted_prob"),
                        labels = c("Actual","Fitted"))

p2 <- ggplot(df2_long, aes(x = value, y = Probability, 
                           color = Type, linetype = Type)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("Actual" = "brown", "Fitted" = "blue")) +
  scale_linetype_manual(values = c("Actual" = "solid", "Fitted" = "dashed")) +
  labs(x = "Value", y = "Probability") +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(5, 5, 5, 5))

## ---------- Combine Side by Side ----------
final_plot <- p1 | p2
final_plot
