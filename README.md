# Computational-Econometrics-R

A light tutotial page on various aspects related to Computational Econometrics issues. 

In many of the following computational and estimation problems due to possibly large computational time when using Statistical Programming Software (such as R) for implementation purposes, it is necessary to use parallel programming techniques to reduce execution time.  

# [A]. Optimization Techniques and Dimensionality Reduction

## [A1]. Portfolio Optimization Techniques 

## References

- Katsouris, C. (2021). Optimal Portfolio Choice and Stock Centrality for Tail Risk Events. [arXiv preprint:2112.12031](https://arxiv.org/abs/2112.12031)
- Brandt, M. W., & Santa‐Clara, P. (2006). Dynamic portfolio selection by augmenting the asset space. The journal of Finance, 61(5), 2187-2217.
- Markowitz, H.M. 1952. Portfolio selection. Journal of Finance 7(1), March, 77–91.
- Markowitz, H.M. 1956. The optimization of a quadratic function subject to linear constraints. Naval Research Logistics Quarterly 3, 111–33.
- Markowitz, H.M. 1959. Portfolio Selection: Efficient Diversification of Investments. New Haven: Yale University Press. Reprinted, New York: John Wiley and Sons, 1970.

## [A2]. Principal Component Analysis 



# [B]. Iterative Simulation and Estimation Methodologies

## [B1.] Minimum Distance Estimation Method

A particular class of econometric models has parameters which are time-varying. For example, consider the case when evaluating the portfolio performance within an out-of-sample framework in which case one has to estimate the covariance matrix in each moving window. Then, a suitable econometric model to capture these time-varying moment conditions is given by 
$$y_t = g_t ( b_0 ) + u_t, \ \ \ \ \ \ \ \ \  t = 1,2,...,T$$

where the response variable yt is an ( m x 1) vector of observable random functions of discrete time and gt is a vector of known functions which depend on a (p x 1) vector of unknown parameters whose true value is denoted by $b0$. In general, gt if a function of a number of exogenous variables as well as b0 so that the function is time dependent. The last component of the model is the vector ut of additive stochastic disturbances. In this case, an appropriate estimation methodology is to employ the iterated minimum distance estimator (MDE). Furthermore, it has been proved by  Phillips (1976) that there is an equivalent relation between the MDE estimator and the quasi-maximum likelihood (QML) estimator when the model disturbances are serially independent.  

Therefore, we consider any vector g ( T ) in which minimizes the following quadratic form 
$$ R_T( b ) = T^{-1} \sum_{t=1}^T \left( y_T - g_t (b) \right)^{\prime} S \left( y_T - g_t (b) \right),$$

given the observations ( yt : t = 1,...,T ) and some positive definite matrix S, the above expression is called a MDE of b0. Then, concentrating the likelihood function with respect to b, we find that bT (S) minimizes the following expression
$$\text{log det} \ \left( T^{-1} \sum_{t=1}^T \left( y_T - g_t (b) \right)^{\prime} \left( y_T - g_t (b) \right) \right).$$


## References

- Wolfowitz, J. (1957). The minimum distance method. The Annals of Mathematical Statistics, 75-88.
- Phillips, P. C. B. (1976). The iterated minimum distance estimator and the quasi-maximum likelihood estimator. Econometrica: Journal of the Econometric Society, 449-460.


## References

- Andrews, D. W., & Lu, B. (2001). Consistent model and moment selection procedures for GMM estimation with application to dynamic panel data models. Journal of econometrics, 101(1), 123-164.
- Hansen, L. (1982): “Large Sample Properties of Generalized Method of Moments Estimators,” Econometrica
- Wu, C. J. (1983). On the convergence properties of the EM algorithm. The Annals of statistics, 95-103.


# Further Reading

[1] Abadir, K. M., & Magnus, J. R. (2005). Matrix algebra (Vol. 1). Cambridge University Press.

[2] Hayashi, F. (2011). Econometrics. Princeton University Press.

[3] White, H. (1996). Estimation, inference and specification analysis (No. 22). Cambridge university press.

[4] Zaslavski, A. J. (2010). Optimization on metric and normed spaces (Vol. 44). Springer Science & Business Media.
