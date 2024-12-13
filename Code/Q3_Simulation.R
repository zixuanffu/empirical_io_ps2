rm(list = ls())
pacman::p_load(evd, data.table, xtable)
load("Data/Q2_NumericalSol.RData")

# simulate state transition
# we pick initial state as (0,0,2)
n <- 10000
idx <- 6
s_init <- s[idx, ]
epsilon_0 <- rgumbel(n)
epsilon_1 <- rgumbel(n)

# start the loop
beta <- 0.99
M <- list(M_0, M_1)
s_next <- s_init

s_list <- as.data.frame(NULL, col.names = c("i", "c", "p"))
x_list <- c()
idx_list <- c()

for (i in 1:n) {
    s_list <- rbind(s_list, s_next)
    idx_list <- c(idx_list, idx)
    x <- ((u(as.matrix(s_next), 1) + epsilon_1[i] + beta * (M[[2]][idx, ] %*% V_ss)) > (u(as.matrix(s_next), 0) + epsilon_0[i] + beta * (M[[1]][idx, ] %*% V_ss)))
    pr <- M[[x + 1]][idx, ]
    idx <- sample(1:20, 1, prob = M[[x + 1]][idx, ])
    s_next <- s[idx, ]
    x_list <- c(x_list, x)
}

saveRDS(s_list, "Data/Q3_s_list.RDS")
saveRDS(x_list, "Data/Q3_x_list.RDS")

dt <- cbind(s_list, x_list)
setnames(dt, c("i", "c", "p", "x"))
setDT(dt)
saveRDS(dt, "Data/Q3_Simulation.RDS")

# compute summary statistics from simulation
dt <- readRDS("Data/Q3_Simulation.RDS")
x_freq <- dt[x == 1, .N] / n
prob_x_lowp <- dt[p == 0.5 & x == 1, .N] / dt[p == 0.5, .N]
dur_lowp <- n / dt[p == 0.5, .N]
dur_x <- n / dt[x == 1, .N]

name <- c("Frequency of purchase", "Probability of purchase when sales", "Average duration between sales", "Average duration between purchases")
des <- data.table(name, c(x_freq, prob_x_lowp, dur_lowp, dur_x))
setnames(des, c("statistic", "value"))
print(xtable(des, digits = 3), floating = FALSE, type = "latex", file = "Results/Tables/simulation_des.tex")

# compute summary statistics from numerical solution
# calculate the conditional choice probability
ccp <- function(u_0, u_1, M_0, M_1, V, beta = 0.99) {
    res <- exp(u_1 + beta * M_1 %*% V) / (exp(u_0 + beta * M_0 %*% V) + exp(u_1 + beta * M_1 %*% V))
    return(res)
}

prob_x <- ccp(u_0, u_1, M_0, M_1, V_ss) # conditional choice probability
pi_1 <- M_1 * matrix(rep(prob_x, 20), ncol = 20)
pi_0 <- M_0 * matrix(rep(1 - prob_x, 20), ncol = 20)
pi_s <- pi_1 + pi_0 # the transition matrix

# find the stationary distribution is to find the eigenvector of pi transpose
eigen <- eigen(t(pi_s))
eigen_vector <- eigen$vectors[, 1]
stationary <- eigen_vector / sum(eigen_vector)
x_freq_theo <- sum(prob_x * stationary) # frequency of purchase: 0.552
# note that for the state vector (i, c, p), low price state = state with odd index
stationary_lowp <- stationary[seq(1, 20, 2)] / sum(stationary[seq(1, 20, 2)])
prob_x_lowp_theo <- sum(prob_x[seq(1, 20, 2)] * stationary_lowp) # 0.62
dur_lowp_theo <- 1 / sum(stationary[seq(1, 20, 2)]) # 1.26
dur_x_theo <- 1 / x_freq_theo # 1.81

name <- c("Frequency of purchase", "Probability of purchase when sales", "Average duration between sales", "Average duration between purchases")
des_theo <- data.table(name, c(x_freq_theo, prob_x_lowp_theo, dur_lowp_theo, dur_x_theo))
setnames(des_theo, c("statistic", "value"))
print(xtable(des_theo, digits = 3), floating = FALSE, type = "latex", file = "Results/Tables/theoretical_des.tex")
