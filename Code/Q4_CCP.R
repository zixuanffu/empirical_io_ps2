rm(list = ls())
pacman::p_load(rootSolve, data.table, xtable)

load("Data/Q2_NumericalSol.RData")
dt <- readRDS("Data/Q3_Simulation.RDS")
setDT(s)
s[, id := 1:.N]

# estimate choice probability Pr(x=1|s)
total <- dt[, .N, by = .(i, c, p)]
purchase <- dt[x == 1, .N, by = .(i, c, p)]
total <- merge(total, purchase, by = c("i", "c", "p"))
setnames(total, c("i", "c", "p", "Total", "Purchase"))
total[, Prob_x := Purchase / Total]
total <- merge(total, s, by = c("i", "c", "p"))
hat_ccp <- total[order(id), ]$Prob_x

dt[, rowid := 1:.N]
dt_new <- merge(dt, s, all.x = TRUE)
dt_new <- dt_new[order(rowid), ]
dt_new[, id_next := shift(id, type = "lead")]
dt_new <- dt_new[-10000, ]

# estimate transition probability matrix for x=0
state_history <- dt_new[x == 0, .(rowid, id, id_next)]
transition_history <- state_history[, .N, by = .(id, id_next)]
matrix_count <- dcast(transition_history, id ~ id_next, value.var = "N")
matrix_count[is.na(matrix_count)] <- 0
matrix_count[, row_sum := rowSums(matrix_count)]
matrix_prob <- matrix_count[, lapply(.SD, function(x) x / row_sum), .SDcols = 2:21]
hat_M_0 <- as.matrix(matrix_prob)

saveRDS(hat_M_0, "Data/Q4_hat_M_0.RDS")

# estimate transition probability matrix for x=1
state_history <- dt_new[x == 1, .(rowid, id, id_next)]
transition_history <- state_history[, .N, by = .(id, id_next)]
matrix_count <- dcast(transition_history, id ~ id_next, value.var = "N")
matrix_count[is.na(matrix_count)] <- 0
matrix_count[, row_sum := rowSums(matrix_count)]
matrix_prob <- matrix_count[, lapply(.SD, function(x) x / row_sum), .SDcols = 2:21]
hat_M_1 <- as.matrix(matrix_prob)

saveRDS(hat_M_1, "Data/Q4_hat_M_1.RDS")

# solve for the expected value function $\bar{V}$
hat_exp_v <- function(V, u_0, u_1, M_0, M_1, prob_x, beta = 0.99, gamma = 0.5772157) {
    rhs0 <- (1 - prob_x) * (u_0 - log(1 - prob_x) + beta * (M_0 %*% V))
    rhs1 <- prob_x * (u_1 - log(prob_x) + beta * (M_1 %*% V))
    res <- (V - (gamma + rhs0 + rhs1))
    return(res)
}

hat_exp_v <- function(V, u_0, u_1, M_0, M_1, prob_x, beta = 0.99, gamma = 0.5772157) {
    res <- V - (gamma + u_0 + beta * M_0 %*% V - log(1 - prob_x))
    return(res)
}

# exp_v <- function(V, u_0, u_1, M_0, M_1, beta = 0.99, gamma = 0.5772157) {
#     res <- V - (gamma + log(exp(u_0 + beta * M_0 %*% V) + exp(u_1 + beta * M_1 %*% V)))
#     return(res)
# }

V_init <- rep(1, 20)
solution <- multiroot(hat_exp_v, start = V_init, u_0 = u_0, u_1 = u_1, M_0 = hat_M_0, M_1 = hat_M_1, prob_x = hat_ccp)
hat_V_ss <- solution$root

hat_V_ss
V_ss
V_both <- cbind(V_state, hat_V_ss)
print(xtable(V_both), floating = FALSE, type = "latex", file = "Results/Tables/V_both.tex")


# export the estimated transition matrix
hat_M_0 <- readRDS("Data/Q4_hat_M_0.RDS")
hat_M_1 <- readRDS("Data/Q4_hat_M_1.RDS")

hat_M_0 <- round(hat_M_0, 3)
hat_M_1 <- round(hat_M_1, 3)

M_0_both <- matrix(NA, 20, 20)
M_1_both <- matrix(NA, 20, 20)
for (i in 1:nrow(M_0)) {
    for (j in 1:ncol(M_0)) {
        # Format the values with the second matrix in parentheses
        M_0_both[i, j] <- paste0(M_0[i, j], " (", hat_M_0[i, j], ")")
    }
}
print(xtable(M_0_both), floating = FALSE, type = "latex", file = "Results/Tables/M_0_both.tex")

for (i in 1:nrow(M_1)) {
    for (j in 1:ncol(M_1)) {
        # Format the values with the second matrix in parentheses
        M_1_both[i, j] <- paste0(M_1[i, j], " (", hat_M_1[i, j], ")")
    }
}
print(xtable(M_1_both), floating = FALSE, type = "latex", file = "Results/Tables/M_1_both.tex")
