# Computational-Econometrics-R

A light tutotial page on various aspects related to Computational Econometrics issues. 

In many of the following econometric methods (especially when considering high-dimensional inference problems) which are implemented on Econometric/Statistical Programming Software (such as R), large computational time requires to consider parallel programming techniques in order to reduce execution time.  

Now 'parallel programming' techniques can be employed in various forms depending on the programming environment or the operating system. Firstly, these techniques do not necessarily refer to the econometric estimation method, however many times the executation time can be reduced when the programming algorithm is 'optimized' in terms of the way that operations/functions are executed. Consider for instance, the use of the Bootstrap Monte Carlo Algorithm to obtain critical values for the underline distribution of a test statistic which requires simulation techniques to obtain asymptotic approxiations. Therefore, running R scripts using RStudio on Personal Machines can be speed-up by using build-in packages that perform parallel optimization (see, R package ['parallel'](https://cran.r-project.org/web/views/HighPerformanceComputing.html)). Secondly, if one considers executing R scripts using HPC then 'parallel programming' refers to the optimal execution of those R scripts that use the available computer resources (e.g., computing capabilities) as efficient as possible. See the following bibliography:

- Eugster, M. J., Knaus, J., Porzelius, C., Schmidberger, M., & Vicedo, E. (2011). Hands-on tutorial for parallel computing with R. Computational Statistics, 26(2), 219-239.
- [Parallel Computation in R.](https://bookdown.org/rdpeng/rprogdatascience/parallel-computation.html)

# [A]. Optimization Techniques and Dimensionality Reduction

## [A1]. Portfolio Optimization Techniques 

Consider the traditional minimum-variance portfolio allocation problem

$$\underset{ \boldsymbol{w} }{ \mathsf{arg min} } \ \boldsymbol{w}^{\top} \boldsymbol{\Sigma} \boldsymbol{w} \ \ \ \text{subject to} \ \ \ \boldsymbol{w}^{\top} \boldsymbol{1} = 1.$$

where w = ( w1,..., wp ) represents the weights put on different assets and 1 is a p-dimensional vector with all entries being 1. Then, the optimal vector of weights which has a closed-form solution is given by the following expression: 

$$w^{*} = \frac{ \boldsymbol{\Sigma}^{-1} \boldsymbol{1} }{ \boldsymbol{1}^{\top} \boldsymbol{\Sigma}^{-1} \boldsymbol{1} }.$$

## References

- Katsouris, C. (2021). Optimal Portfolio Choice and Stock Centrality for Tail Risk Events. [arXiv preprint:2112.12031](https://arxiv.org/abs/2112.12031)
- Ledoit, O., & Wolf, M. (2022). Quadratic shrinkage for large covariance matrices. Bernoulli, 28(3), 1519-1547.
- Maillet, B., Tokpavi, S., & Vaucher, B. (2015). Global minimum variance portfolio optimisation under some model risk: A robust regression-based approach. European Journal of Operational Research, 244(1), 289-299.
- Brandt, M. W., & Santa‐Clara, P. (2006). Dynamic portfolio selection by augmenting the asset space. The journal of Finance, 61(5), 2187-2217.
- Markowitz, H.M. 1952. Portfolio selection. Journal of Finance 7(1), March, 77–91.
- Markowitz, H.M. 1956. The optimization of a quadratic function subject to linear constraints. Naval Research Logistics Quarterly 3, 111–33.
- Markowitz, H.M. 1959. Portfolio Selection: Efficient Diversification of Investments. New Haven: Yale University Press. Reprinted, New York: John Wiley and Sons, 1970.
- Yu, K. B. (1991). Recursive updating the eigenvalue decomposition of a covariance matrix. IEEE Transactions on Signal Processing, 39(5), 1136-1145.

### [A1.1.] Testing for Stochastic Dominance

Tests for Stochastic Dominance are usually nonparametric tests with mainy applications when comparing observationally equivalent structures. Under certain regularity conditions stochastic dominance tests can be employed to evaluate the stochastic dominance of portfolios and portfolio returns. In particular, Second-Order Stochastic Dominance can be evaluated within the framework of expected utility maximization. The main idea of such a statistical test is the interest in examining stochastic monotonic relationships for economic and financial relationships. Consider for example, two continuous random variables X and Y both supported on the space [0,1]. Then, the following assumptions should hold for the statistical validiy of these tests.  

1. Investors use an expected utility maximization problem along with the portfolio returns 
2. The asset returns are assumed to be serially independent and identically distributed with 

$$E[ \mathbf{x} ] = \mathbf{\mu} \ \ \ \text{and} \ \ E[ (\mathbf{x} - \mathbf{\mu}).(\mathbf{x} - \mathbf{\mu})^{\top}] = \mathbf{\Sigma}.$$ 

3. Investors are allowed to diversify between assets with a corresponding N-dimensional vector of portfolio weights.

Additionally, in order to be able to evaluate the effectiveness of the portfolio allocation problem and the corresponding portfolio returns induced by a risk matrix such as the financial connectedness matrix proposed by Katsouris (2021), we may also be interested to examine the robustness of the test under the distribution invariance assumption, or under changes in the centrality structure or network topology of assets within the network. For example, the effectiveness of an investment strategy may be better be evaluated under the assumption that there are no linkages or spillover effects between the nodes of the network.


## [A2]. Principal Component Analysis 



## References

- Tipping, M. E., & Bishop, C. M. (1999). Probabilistic principal component analysis. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 61(3), 611-622.
- James, W., & Stein, C. (1961). Estimation with quadratic loss Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability, Volume 1: Contributions to the Theory of Statistics, Berkeley.
- Tikhonov, A. N. (1963). On the solution of ill-posed problems and the method of regularization. In Doklady Akademii Nauk (Vol. 151, No. 3, pp. 501-504). Russian Academy of Sciences.



# [B]. Iterative Simulation and Estimation Methodologies

## [B1.] Minimum Distance Estimation Method

A particular class of econometric models has parameters which are time-varying. For example, consider the case when evaluating the portfolio performance within an out-of-sample framework in which case one has to estimate the covariance matrix in each moving window. Then, a suitable econometric model to capture these time-varying moment conditions is given by 
$$y_t = g_t ( b_0 ) + u_t, \ \ \ \ \ \ \ \ \  t = 1,2,...,n$$

where the response variable yt is an ( m x 1) vector of observable random functions of discrete time and gt is a vector of known functions which depend on a (p x 1) vector of unknown parameters whose true value is denoted by $b0$. In general, gt if a function of a number of exogenous variables as well as b0 so that the function is time dependent. The last component of the model is the vector ut of additive stochastic disturbances. In this case, an appropriate estimation methodology is to employ the iterated minimum distance estimator (MDE). Furthermore, it has been proved by  Phillips (1976) that there is an equivalent relation between the MDE estimator and the quasi-maximum likelihood (QML) estimator when the model disturbances are serially independent.  

Therefore, we consider any vector g ( n ) in which minimizes the following quadratic form 

$$ R_n( b ) = n^{-1} \sum_{t=1}^n ( y_n - g_t (b) )^{\prime} S ( y_n - g_t (b) )$$

given the observations ( yt : t = 1,...,n ) and some positive definite matrix S, the above expression is called a MDE of b0. Then, concentrating the likelihood function with respect to b, we find that bn (S) minimizes the following expression

$$\text{log det} \ \left( n^{-1} \sum_{t=1}^n \left( y_n - g_t (b) \right)^{\prime} \left( y_n - g_t (b) \right) \right).$$

## Example 1

Consider the MDE of an AR(1) model using R. 


```R

## See the R package 'KoulMde'  

n <- 10
p <- 3

## Generate n-by-p design matrix X
X <- matrix(runif(n*p, 0,50), nrow=n, ncol=p)  

## Generate true model parameters
beta <- c(-2, 0.3, 1.5)                      
rho  <- 0.4                                   
eps  <- vector(length=n)
xi  <- rnorm(n, 0,1)                          
                                            
for(i in 1:n)
{
 if(i==1){eps[i] <- xi[i]}
 else{eps[i] <- rho*eps[i-1] + xi[i]}
}
Y <- X\%*\%beta + eps

## Use the default weight matrix
D <- "default"  

## Set initial value for beta
b0 <- solve(t(X)\%*\%X)\%*\%(t(X)\%*\%Y)             

## Define Lebesgue measure
IntMeasure      <- "Lebesgue"                                
MDEResult       <- Koul2StageMde(Y,X, "default", b0, IntMeasure, 1, IntMeasure, TuningConst = 1.345)
MDE1stageResult <- MDEResult[[1]]
MDE2stageResult <- MDEResult[[2]]

beta1     <- MDE1stageResult$betahat1stage
residual1 <- MDE1stageResult$residual1stage
rho1      <- MDE1stageResult$rhohat1stage

beta2     <- MDE2stageResult$betahat2stage
residual2 <- MDE1stageResult$residual2stage
rho2      <- MDE2stageResult$rhohat2stage


```
## References

- Wolfowitz, J. (1957). The minimum distance method. The Annals of Mathematical Statistics, 75-88.
- Phillips, P. C. B. (1976). The iterated minimum distance estimator and the quasi-maximum likelihood estimator. Econometrica: Journal of the Econometric Society, 449-460.
- Koul, H. L. (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. The Annals of Statistics, 1194-1213.


## [B2.] Generalized Method of Moments (GMM)

An alternative estimation technique to MDE and MLE is the Method of Moments, where the main idea is the mapping of sample moments to population moments. Furthermore, in the case of Generalized Method of Moments (GMM) there are more moments than parameters and thus econometric identification can be achieved. In addition, GMM can be employed regardless of any distributional assumptions which makes it an attractive alternative to MLE especially when Gaussianity assumption is not required.     

Define with 

$$g_n(  X, \theta ) = \frac{1}{n}  \sum_{i=1}^n g( X_i; \theta ) \ \ \ \text{such that} \ \ \mathbb{E} [ g_n(  X, \theta_0 ) ] = 0.$$

Then, the corresponding GMM estimator is obtained based on the following optimization function 

$$\hat{\theta}_n = \underset{  \theta \in \Theta  }{ \mathsf{argmin}  }  \ \ [ g_n(  X, \theta ) ]^{\prime} W [ g_n(  X, \theta ) ] $$

In other words, the statistical geometry of the GMM estimator minimizes the distance of sample counterparts to popultation parameters based on the underline orthogonality conditions and the distance matrix W.  


## Example 2


## References

- Andrews, D. W., & Lu, B. (2001). Consistent model and moment selection procedures for GMM estimation with application to dynamic panel data models. Journal of econometrics, 101(1), 123-164.
- Chesher, A., & Smith, R. J. (1997). Likelihood ratio specification tests. Econometrica: Journal of the Econometric Society, 627-646.
- Hansen, L. (1982): “Large Sample Properties of Generalized Method of Moments Estimators,” Econometrica
- Wu, C. J. (1983). On the convergence properties of the EM algorithm. The Annals of statistics, 95-103.


# Further Reading

[1] Abadir, K. M., & Magnus, J. R. (2005). Matrix algebra (Vol. 1). Cambridge University Press.

[2] Hayashi, F. (2011). Econometrics. Princeton University Press.

[3] White, H. (1996). Estimation, inference and specification analysis (No. 22). Cambridge University Press.

[4] Boyd, S., Boyd, S. P., & Vandenberghe, L. (2004). Convex optimization. Cambridge University Press.

[5] Zaslavski, A. J. (2010). Optimization on metric and normed spaces (Vol. 44). Springer Science & Business Media.

# Disclaimer

The author (Christis G. Katsouris) declares no conflicts of interest.

The proposed Course Syllabus is currently under development and has not been officially undergone quality checks. All rights reserved.
