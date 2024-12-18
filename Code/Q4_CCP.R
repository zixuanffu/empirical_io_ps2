rm(list = ls())
pacman::p_load(rootSolve, data.table, xtable)

load("Data/Q2_NumericalSol.RData")
dt <- readRDS("Data/Q3_Simulation.RDS")
setDT(s)
s[, id := 1:.N]

# method 1: estimate choice probability Pr(x=1|s) by frequency estimator
total <- dt[, .(Total = .N, Purchase = sum(x)), by = .(i, c, p)]
total[, Prob_x := Purchase / Total]
total <- merge(s, total, all.x = TRUE) # we need ccp for all states, what should we do if some states are not visited?
total[is.na(Prob_x)]$id # 62 66 65 68 67
total[Prob_x == 0]$id # 60 61
total[Prob_x == 1]$id # 2 46 50 54 53 56 58 57 64
total[is.na(Prob_x), Prob_x := 0]
hat_ccp <- total[order(id), Prob_x]


# method 2: estimate by kernel density estimator with Gaussian kernel

# bandwidth selection

bandwidth <- 0.5

# gaussian kernel
gaussian_kernel <- function(d, bw) {
    exp(-0.5 * (d / bw)^2) / (bw * sqrt(2 * pi))
}

# smoothing
for (row in 1:nrow(s)) {
    target_state <- s[row, .(i, c, p)] # Current state

    distances <- sqrt((dt$i - target_state$i)^2 +
        (dt$c - target_state$c)^2 +
        (dt$p - target_state$p)^2) # Euclidean distance

    weights <- gaussian_kernel(distances, bandwidth) # Apply kernel
    numerator <- sum(weights * dt$x) # Weighted sum of choices
    denominator <- sum(weights) # Sum of weights

    # Estimate smoothed probability
    s[row, Prob_x := ifelse(denominator > 0, numerator / denominator, 0)]
}

hat_ccp <- s[order(id), Prob_x]

# ccp estimator
hat_exp_v <- function(V, u_0, u_1, M_0, M_1, prob_x, beta = 0.99, gamma = 0.5772157) {
    res <- V - (gamma + (1 - prob_x) * (-log(1 - prob_x) + u_0 + beta * M_0 %*% V) + prob_x * (-log(prob_x) + u_1 + beta * M_1 %*% V))
    return(res)
}

V_init <- rep(1, nrow(s))
solution <- multiroot(hat_exp_v, start = V_init, u_0 = u_0, u_1 = u_1, M_0 = M_0, M_1 = M_1, prob_x = hat_ccp)
hat_V_ss <- solution$root

hat_V_ss
V_ss
sd(hat_V_ss)
sd(V_ss)

V_both <- cbind(V_state, hat_V_ss)
colnames(V_both) <- c("Inventory", "Consumer purchase", "Price", "$\\bar{V}$", "$\\bar{V}_{ccp}$")
print(xtable(V_both), floating = FALSE, type = "latex", file = "Results/Tables/V_both.tex", sanitize.colnames.function = identity)
