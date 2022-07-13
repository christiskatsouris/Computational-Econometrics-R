# Computational-Econometrics-R

A light tutotial page on various aspects related to Computational Econometrics issues. 

In many of the following econometric methods (especially when considering high-dimensional inference problems) which are implemented on Econometric/Statistical Programming Software (such as R), large computational time requires to consider parallel programming techniques in order to reduce execution time.  

Now 'parallel programming' techniques can be employed in various forms depending on the programming environment or the operating system. Firstly, these techniques do not necessarily refer to the econometric estimation method, however many times the executation time can be reduced when the programming algorithm is 'optimized' in terms of the way that operations/functions are executed. Consider for instance, the use of the Bootstrap Monte Carlo Algorithm to obtain critical values for the underline distribution of a test statistic which requires simulation techniques to obtain asymptotic approxiations. Therefore, running R scripts using RStudio on Personal Machines can be speed-up by using build-in packages that perform parallel optimization (see, R package ['parallel'](https://cran.r-project.org/web/views/HighPerformanceComputing.html)). Secondly, if one considers executing R scripts using HPC then 'parallel programming' refers to the optimal execution of those R scripts that use the available computer resources (e.g., computing capabilities) as efficient as possible. See the following bibliography:

- Eugster, M. J., Knaus, J., Porzelius, C., Schmidberger, M., & Vicedo, E. (2011). Hands-on tutorial for parallel computing with R. Computational Statistics, 26(2), 219-239.
- [Parallel Computation in R.](https://bookdown.org/rdpeng/rprogdatascience/parallel-computation.html)

```R


```

# [A]. Optimization Techniques and Dimensionality Reduction

## [A1]. Portfolio Optimization Techniques 

Consider rN to be the N-dimensional vector of stock returns at time t. In other words, this is considered to be a strictly stationary time series of an (N x 1) observed vector of returns (e.g., log returns) with population moments given by 

$$ \boldsymbol{\mu} = \mathbb{E}[ \boldsymbol{r}_t ], \ \ \text{where} \ \boldsymbol{\mu} \in \mathbb{R}^{N \times 1}.$$

representing the vector of expected returns for a given time period, e.g., a historical period or a full sample period $t=1,...,T$. 

$$\boldsymbol{\Sigma}= \mathsf{Cov}[ \boldsymbol{r}_t ] = \mathbb{E}[ \boldsymbol{r}_t \boldsymbol{r}_t^{\top} ] - \boldsymbol{\mu} \boldsymbol{\mu}^{\top}.$$
representing the (unconditional) covariance matrix (often called volatility matrix) of stock returns for the time period $t=1,...,T$. 

Moreover, let R_t+1 denote the continuous random variable which represents the portfolio returns at time t+1 which is given by
$$R_{t+1} = \sum_{i=1}^N w_{i,t} r_{i,t+1} \equiv \boldsymbol{w}^{\top} \boldsymbol{r}_{t+1}.$$

Then, the traditional minimum-variance portfolio allocation problem

$$\underset{ \boldsymbol{w} \in \mathbb{R}^N }{ \mathsf{arg min} } \ \boldsymbol{w}^{\top} \boldsymbol{\Sigma} \boldsymbol{w} \ \ \ \text{subject to} \ \ \ \boldsymbol{w}^{\top} \boldsymbol{1} = 1.$$

where w = ( w1,..., wN ) represents the weights put on different assets and 1 is a N-dimensional vector with all entries being 1. Then, the optimal vector of weights which has a closed-form solution is given by the following expression: 

$$w^* = \frac{ \boldsymbol{\Sigma}^{-1} \boldsymbol{1} }{ \boldsymbol{1}^{\top} \boldsymbol{\Sigma}^{-1} \boldsymbol{1} }.$$

Therefore, for high-dimensional portfolio optimization problems, a challenging problem is the robust estimation of the inverse covariance matrix (precision matrix). In practise since the sample covariance matrix is used as a statistic of the corresponding population covariance matrix, then the resulting 'plug-in' portfolio vector has been widely adopted. However a related research question of interest is: How well does such a portfolio perform? A first step in this direction is to consider evaluating the portfolio performance - especially when comparing different estimates, using for instance stochastic dominance tests and/or via other portfolio performance measures. A second step in tackling this aspect is to investigate the asymptotic properties of the portfolio using different estimators for the covariance matrix by deriving appropriate error bounds in relation to different optimization constraints. 

Furthermore, the mean-variance portfolio optimization problem is defined as below

$$ \boldsymbol{w}^* =  \underset{ \boldsymbol{w} \in \mathbb{R}^N  }{ \mathsf{arg min} } \ \boldsymbol{w}^{\top} \boldsymbol{\Sigma} \boldsymbol{w}  \ \ \ \text{subject to} \ \ \ \boldsymbol{w}^{\top} \boldsymbol{\mu} = \boldsymbol{\mu}_0.$$

The closed-form solution of the above investement strategy is given by

$$\boldsymbol{w}^* = \frac{ \boldsymbol{ \mu}_0 }{ \boldsymbol{\mu}^{\top} \boldsymbol{\Sigma}^{-1} \boldsymbol{\mu} } \boldsymbol{\Sigma}^{-1} \boldsymbol{\mu}.$$

Similarly we can also consider Expected Utility Strategies, in which we incorporate the expected utility of an investor and thus the degree of risk aversion within the objective function. 

$$\underset{ \boldsymbol{w} \in \mathbb{R}^N  }{ \mathsf{arg max} } \ U( R_{t+1} ) =  \boldsymbol{w}^{\top} R_{t+1} -  \boldsymbol{w}^{\top} \boldsymbol{\Sigma}_{t+1|t} \boldsymbol{w}$$ 

with optimization constraints given by 

$$\boldsymbol{w}_t^{\top} \boldsymbol{1}= 1 \ \ \ \text{and} \ \ \  0 \leq \boldsymbol{w}_t \leq 1$$

## Example 1

To solve the optimal portfolio problem we solve the Lagrangian function with only constraint the full investment condition. Construct an R coding procedure that obtaints the closed-form soluton based on the Markowitz portfolio optimization problem. Hint: To avoid the presence of ill-conditioned precision matrices especially when considering high-dimensional portfolios, use the pseudo-inverse methodology. 

```R

optimal_weights_function_foc <- function( N = N, Sigma = Sigma )
{#begin of function
  
  #Assign the input values of the function
  N     <- N
  Sigma <- Sigma
  
  pinv <- function(A, eps=1e-8)
  {
    L <- svd(A)
    d <- L$d
    i <- abs(d) > eps
    d[i] <- 1/d[i]
    L$v %*% diag(d, nrow=length(d)) %*% t(L$u)
  }
  
  ones <- matrix(1, nrow = N, ncol = 1)
  sigma_inverse <- pinv(Sigma)
  numerator     <- sigma_inverse%*%ones
  denominator   <- ( t(ones) %*% (sigma_inverse) %*% (ones) )
  denominator   <- as.numeric(denominator)
  
  optimal.weights <- (1/denominator)*numerator
  optimal.weights <- as.matrix(optimal.weights)
  
  # Normalize Weights  
  optimal.weights.norm = optimal.weights/sum(optimal.weights)
  return(optimal.weights)
  
}#end of function

```

In practise, the additional condition of having the weight vector being between zero and one is necessary when constructing the optimization problem to avoid infeasible solutions especially in the case where the optimal portfolio problem is constructed based on small number of observations. To implement the aformentioned procedure in R, we can employ a specific type of optimization algorithms called ['Genetic Algorithm'](https://rpubs.com/Argaadya/550805) (GA) which employs an iterative procedure to converge to an approximate solution given a prespecified number of iterations.  

```R

### R Package for the optimization
install.packages("GA")
install.packages("Matrix")
install.packages("matlib")

library(GA)
library(Matrix)
library(matlib)

optimal_weights_function_GA <- function( N = N, Sigma = Sigma )
{# begin of function
  
  set.seed(1234)
  # Assign the input values of the function
  N     <- N
  Sigma <- Sigma
  
  # normalised weights
  weights <- function(w) 
  { 
    drop(w/sum(w)) 
  }
  
  # Estimation of the variance of the portfolio
  VarPortfolio <- function( w, S=Sigma ) 
  {
    S <- S
    w <- weights(w)
    drop(t(w) %*% S %*% w) # Objective function
  }
  
  # fitness function
  fitness <- function(w) 
  {# begin of function
    constraint = function(w) 
    {# begin of function
      boundary_constr = (sum(w)-1)^2   # "sum x = 1" constraint
      
      for (i in 1:length(w))
      {
        boundary_constr = boundary_constr + 
          max(c(0,w[i]-1))^2 +  # "x <= 1" constraint
          max(c(0,-w[i]))^2     # "x >= 0" constraint
      }
      
     return (boundary_constr)
    }# end of function
    
    penalty <- constraint(w)
    -(VarPortfolio(w) + penalty)
  }# end of function
  
  # Genetic Algorithm Estimation of optimal weights
  GA <- ga(type = "real-valued", fitness = fitness,
           lower = rep(0, Nr_C), upper = rep(1, Nr_C), 
           maxiter = 1000, run = 200, optim = TRUE)
  
  optimal.weights <- weights(GA@solution)
  return( optimal.weights )
  
}# end of function

```

## References

- Katsouris, C. (2021). Optimal Portfolio Choice and Stock Centrality for Tail Risk Events. [arXiv preprint:2112.12031](https://arxiv.org/abs/2112.12031)
- Fan, J., Liao, Y., & Shi, X. (2015). Risks of large portfolios. Journal of Econometrics, 186(2), 367-387.
- Maillet, B., Tokpavi, S., & Vaucher, B. (2015). Global minimum variance portfolio optimisation under some model risk: A robust regression-based approach. European Journal of Operational Research, 244(1), 289-299.
- Brandt, M. W., & Santa‐Clara, P. (2006). Dynamic portfolio selection by augmenting the asset space. The journal of Finance, 61(5), 2187-2217.
- Markowitz, H.M. 1952. Portfolio selection. Journal of Finance 7(1), March, 77–91.
- Markowitz, H.M. 1956. The optimization of a quadratic function subject to linear constraints. Naval Research Logistics Quarterly 3, 111–33.
- Markowitz, H.M. 1959. Portfolio Selection: Efficient Diversification of Investments. New Haven: Yale University Press. Reprinted, New York: John Wiley and Sons, 1970.
- Senneret, M., Malevergne, Y., Abry, P., Perrin, G., & Jaffres, L. (2016). Covariance versus precision matrix estimation for efficient asset allocation. IEEE Journal of Selected Topics in Signal Processing, 10(6), 982-993.
- Yu, K. B. (1991). Recursive updating the eigenvalue decomposition of a covariance matrix. IEEE Transactions on Signal Processing, 39(5), 1136-1145.

## Further Reading

- Antell, J., & Vaihekoski, M. (2019). Expected and realized returns in conditional asset pricing models: A new testing approach. Journal of Empirical Finance, 52, 220-236.
- Avella-Medina, M., Battey, H. S., Fan, J., & Li, Q. (2018). Robust estimation of high-dimensional covariance and precision matrices. Biometrika, 105(2), 271-284.
- Ledoit, O., & Wolf, M. (2022). Quadratic shrinkage for large covariance matrices. Bernoulli, 28(3), 1519-1547.
- Tyler, D. E., & Yi, M. (2020). Lassoing eigenvalues. Biometrika, 107(2), 397-414.
- Zhao, Z., Ledoit, O., & Jiang, H. (2020). Risk reduction and efficiency increase in large portfolios: leverage and shrinkage. University of Zurich, Department of Economics, Working Paper, (328).



### [A1.1.] Testing for Stochastic Dominance

Tests for Stochastic Dominance are usually nonparametric tests with mainy applications when comparing observationally equivalent structures. Under certain regularity conditions stochastic dominance tests can be employed to evaluate the stochastic dominance of portfolios and portfolio returns. In particular, Second-Order Stochastic Dominance can be evaluated within the framework of expected utility maximization. The main idea of such a statistical test is the interest in examining stochastic monotonic relationships for economic and financial relationships. Consider for example, two continuous random variables X and Y both supported on the space [0,1]. Then, the following assumptions should hold for the statistical validiy of these tests.  

(i).   Investors use an expected utility maximization problem along with the portfolio returns. 
(ii).  The asset returns are assumed to be serially independent and identically distributed with 

$$E[ \mathbf{x} ] = \mathbf{\mu} \ \ \ \text{and} \ \ E[ (\mathbf{x} - \mathbf{\mu}).(\mathbf{x} - \mathbf{\mu})^{\top}] = \mathbf{\Sigma}.$$ 

(iii). Investors are allowed to diversify between assets with a corresponding N-dimensional vector of portfolio weights.

Additionally, in order to be able to evaluate the effectiveness of the portfolio allocation problem and the corresponding portfolio returns induced by a risk matrix such as the financial connectedness matrix proposed by [Katsouris (2021)](https://arxiv.org/abs/2112.12031), we may also be interested to examine the robustness of the test under the distribution invariance assumption, or under changes in the centrality structure or network topology of assets within the network. For example, the effectiveness of an investment strategy is feasible to be evaluated under the assumption that there are no linkages or spillover effects between the nodes of the network.

### References
- Chou, P. H., & Zhou, G. (2006). Using Bootstrap to Test Portfolio Efficiency. Annals of Economics & Finance, 7(2).

## [A2]. Principal Component Analysis 



## References

- Tipping, M. E., & Bishop, C. M. (1999). Probabilistic principal component analysis. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 61(3), 611-622.
- James, W., & Stein, C. (1961). Estimation with quadratic loss Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability, Volume 1: Contributions to the Theory of Statistics, Berkeley.
- Tikhonov, A. N. (1963). On the solution of ill-posed problems and the method of regularization. In Doklady Akademii Nauk (Vol. 151, No. 3, pp. 501-504). Russian Academy of Sciences.


## Task 1

Following the following steps construct a small Monte Carlo simulation study where the portfolio performance is evaluated using different constraints (thus changing the investment strategy).

- Step 1: Choose an appropriate optimization methodology and set-up the constraint optimization problem with the necessary conditions. 
- Step 2: Compare the convergence properties of the optimization algorithm and the performance of the optimal portfolio using variance-covariance matrix versus using a GARCH-based covariance matrix. 
- Step 3: Simulate from the Multivariate normal distribution the elements of the covariance matrix.
- Step 4: Repeat the procedure B times for the Monte Carlo step and obtain the empirical variance and MSE in order to compare the different optimization methods.

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
- Andrews, D. W. (1997). A stopping rule for the computation of generalized method of moments estimators. Econometrica: Journal of the Econometric Society, 913-931.
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
