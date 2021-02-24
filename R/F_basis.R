#-----------------------------------------------------

tpower <- function(x, t, p) {
    # Truncated p-th power function
    (x - t)^p * (x > t)
}

bbase <- function(x, xdat = NULL, xl = min(x), xr = max(x), deg = 3, dx = 0.1) {
    # Construct B-spline basis
    if (is.null(xdat)) 
        xdat = x
    
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
    B <- (-1)^(deg + 1) * P %*% t(D)
    return(list(B.ik = B, knots.k = knots))
}

