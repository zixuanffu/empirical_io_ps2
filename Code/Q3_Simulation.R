rm(list = ls())
pacman::p_load(evd, data.table, xtable)
load("Data/Q2_NumericalSol.RData")

num_state <- nrow(s)

# simulate state transition
# we pick initial state as (0,0,0.5)
n <- 10000
idx <- 2
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

# simulation
# method 1 (the big transition matrix)
# for (i in 1:n) {
#     s_list <- rbind(s_list, s_next)
#     idx_list <- c(idx_list, idx)
#     x <- ((u(as.matrix(s_next), 1) + epsilon_1[i] + beta * (M[[2]][idx, ] %*% V_ss)) > (u(as.matrix(s_next), 0) + epsilon_0[i] + beta * (M[[1]][idx, ] %*% V_ss)))
#     pr <- M[[x + 1]][idx, ]
#     idx <- sample(1:num_state, 1, prob = M[[x + 1]][idx, ])
#     s_next <- s[idx, ]
#     x_list <- c(x_list, x)
# }


# method 2: individual transition matrix! (m_i0, m_i1, m_c, m_p)
for (i in 1:n) {
    s_list <- rbind(s_list, s_next)
    idx_list <- c(idx_list, idx)
    x <- ((u(as.matrix(s_next), 1) + epsilon_1[i] + beta * (M[[2]][idx, ] %*% V_ss)) > (u(as.matrix(s_next), 0) + epsilon_0[i] + beta * (M[[1]][idx, ] %*% V_ss)))
    # the only difference is the sampling of the next state
    p_cur <- as.numeric(s_next[3])
    i_cur <- as.numeric(s_next[1])

    p_next <- sample(c(2, 0.5), 1, prob = m_p[(p_cur == 0.5) + 1, ])
    c_next <- sample(c(0, 0.25), 1, prob = c(0.5, 0.5))
    if (x) {
        i_next <- sample(seq(0, 4, 0.25), 1, prob = m_i1[i_cur[1] * 4 + 1, ])
    } else {
        i_next <- sample(seq(0, 4, 0.25), 1, prob = m_i0[i_cur[1] * 4 + 1, ])
    }

    idx <- which(s[, 1] == i_next & s[, 2] == c_next & s[, 3] == p_next)
    s_next <- s[idx, ]
    x_list <- c(x_list, x)
}

dim(unique(s_list)) # check the number of unique states that we managed to visit. 65
setDT(s)
s[, id := 1:.N] # add id to the state space
s_unique <- merge(unique(s_list), s, by = c("i", "c", "p"))
setdiff(s$id, s_unique$id) # check the id of the states that we did not visit:  62 65 66 67 68

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
pi_1 <- M_1 * matrix(rep(prob_x, num_state), ncol = num_state)
pi_0 <- M_0 * matrix(rep(1 - prob_x, num_state), ncol = num_state)
pi_s <- pi_1 + pi_0 # the transition matrix

# find the stationary distribution is to find the eigenvector of pi transpose
eigen <- eigen(t(pi_s))
eigen_vector <- Re(eigen$vectors[, 1])
stationary <- eigen_vector / sum(eigen_vector)
s$ss_dist <- stationary # for some states, the long term probability is zero
x_freq_theo <- sum(prob_x * stationary) # frequency of purchase: 0.499
# note that for the state vector (i, c, p), low price state = state with even index
stationary_lowp <- stationary[seq(2, num_state, 2)] / sum(stationary[seq(2, num_state, 2)])
prob_x_lowp_theo <- sum(prob_x[seq(2, num_state, 2)] * stationary_lowp) # 0.762
dur_lowp_theo <- 1 / sum(stationary[seq(2, num_state, 2)]) # 4.8
dur_x_theo <- 1 / x_freq_theo # 2

name <- c("Frequency of purchase", "Probability of purchase when sales", "Average duration between sales", "Average duration between purchases")
des_theo <- data.table(name, c(x_freq_theo, prob_x_lowp_theo, dur_lowp_theo, dur_x_theo))
setnames(des_theo, c("statistic", "value"))
print(xtable(des_theo, digits = 3), floating = FALSE, type = "latex", file = "Results/Tables/theoretical_des.tex")
