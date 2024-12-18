rm(list = ls())
pacman::p_load(rootSolve, xtable)

# discrete state space
i <- seq(0, 4, 0.25)
c <- c(0, 0.25) # 0 or 0.25
p <- c(2, 0.5)
s <- expand.grid(p = p, c = c, i = i)
s <- s[, c(3, 2, 1)]

# transition matrix for i, c, p separately
m_c <- matrix(0.5, 2, 2)
m_p <- matrix(c(0.75, 0.25, 0.95, 0.05), 2, 2, byrow = TRUE)
m_i0 <- diag(0.5, length(i))
m_i0[row(m_i0) - col(m_i0) == 1] <- 0.5
m_i1 <- t(m_i0)
m_i0[1, 1] <- 1
m_i1[length(i), length(i)] <- 1

# overall transition matrix using Kronecker product
M_0 <- kronecker(m_i0, kronecker(m_c, m_p))
M_1 <- kronecker(m_i1, kronecker(m_c, m_p))
# check if the sum of each row is 1
View(rowSums(M_0))
View(rowSums(M_1))
# utility function
u <- function(s, x, lambda = 3, alpha = 2) {
    i <- s[1]
    c <- s[2]
    p <- s[3]
    res <- (-lambda * (c > 0) * (i == 0) + alpha * c - x * p)
    return(as.matrix(res))
}
u_0 <- u(s, 0)
u_1 <- u(s, 1)

# expected value function

# exp_v <- function(V, u_0, u_1, M_0, M_1, beta = 0.99, gamma = 0.5772157) {
#     res <- V - (gamma + u_0 + beta * M_0 %*% V + log(1 + exp((u_1 + beta * M_1 %*% V) - (u_0 + beta * M_0 %*% V))))
#     return(res)
# }

exp_v <- function(V, u_0, u_1, M_0, M_1, beta = 0.99, gamma = 0.5772157) {
    res <- V - (gamma + log(exp(u_1 + beta * M_1 %*% V) + exp(u_0 + beta * M_0 %*% V)))
    return(res)
}

# solve for the $\bar{V}$
V_init <- rep(1, nrow(s))
solution <- multiroot(exp_v, start = V_init, u_0 = u_0, u_1 = u_1, M_0 = M_0, M_1 = M_1)
V_ss <- solution$root
V_ss
V_state <- cbind(s, V_ss)

# save the results

save.image(file = "Data/Q2_NumericalSol.RData")
colnames(V_state) <- c("Inventory", "Consumer purchase", "Price", "Expected value function")
print(xtable(V_state), floating = FALSE, type = "latex", file = "Results/Tables/V_state.tex")
