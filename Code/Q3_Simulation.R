rm(list = ls())
pacman::p_load(evd, data.table)
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

# describe the simulation
dt <- readRDS("Data/Q3_Simulation.RDS")
x_freq <- dt[x == 1, .N] / n
prob_x_lowp <- dt[p == 0.5 & x == 1, .N] / dt[p == 0.5, .N]
dur_lowp <- n / dt[p == 0.5, .N]
dur_x <- n / dt[x == 1, .N]

name <- c("Frequency of purchase", "Probability of purchase when sales", "Average duration between sales", "Average duration between purchases")
des <- data.table(name, c(x_freq, prob_x_lowp, dur_lowp, dur_x))
setnames(des, c("statistic", "value"))
print(xtable(des), floating = FALSE, type = "latex", file = "Results/Tables/simulation_des.tex")
