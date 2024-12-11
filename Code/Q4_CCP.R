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

hat_exp_v <- function(V, u_0, u_1, M_0, M_1, prob_x, beta = 0.99, gamma = 0.5772157) {
    res <- V - (gamma + (1 - prob_x) * (u_0 + beta * M_0 %*% V) + prob_x * (u_1 + beta * M_1 %*% V))
    return(res)
}

V_init <- rep(1, 20)
solution <- multiroot(hat_exp_v, start = V_init, u_0 = u_0, u_1 = u_1, M_0 = M_0, M_1 = M_1, prob_x = hat_ccp)
hat_V_ss <- solution$root

hat_V_ss
V_ss

V_both <- cbind(V_state, hat_V_ss)
colnames(V_both) <- c("Inventory", "Consumer purchase", "Price", "$\\bar{V}$", "$\\bar{V}_{ccp}$")
print(xtable(V_both), floating = FALSE, type = "latex", file = "Results/Tables/V_both.tex", sanitize.colnames.function = identity)
