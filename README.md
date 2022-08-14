# Computational-Econometrics-R

Teaching page on various aspects related to Computational Econometrics issues using the Statistical Software R (Drafted July 2022). 

### Course Overview: 

The main philosophy with this course is to combine traditional statistical learning methods such as the portfolio optimization problem with computational and numerical methods such as the Generalized Method of Momements which are useful for econometric modelling of large datasets. Emphasis will be given to some examples of  time series regression applications that are useful in related computational problems. The module will introduce state-of-the-art techniques and programming software capabilities for each theme covered in the Syllabus. 

> Further extensions can be also considered within a high dimensional environment although these aspects correspond to a course on 'Advanced Computational Econometrics' (to be developed in the future).

## R Parallel Programming Techniques  

In many of the following econometric methods (especially when considering high-dimensional inference problems or when implementing resampling methodologies such as permutation or bootstrap testing) with Econometric/Statistical Programming Software (such as R, Matlab or Stata) in order to reduce the execution time requires to apply parallel programming techniques and/or distributing tasks using HPC.   

Now 'parallel programming' techniques can be employed diffently depending on the programming environment as well as on the structure of the algorithm under consideration. Firstly, these techniques do not necessarily refer to the econometric estimation method, however many times the executation time can be reduced when the programming algorithm is 'optimized' in terms of the way that operations/functions are executed. Consider for instance, the use of the Bootstrap Monte Carlo Algorithm to obtain critical values for the underline distribution of a test statistic which requires simulation techniques to obtain asymptotic approxiations. Therefore, running R scripts using RStudio on Personal Machines can be speed-up by using build-in packages that perform parallel optimization (see, R package ['parallel'](https://cran.r-project.org/web/views/HighPerformanceComputing.html)). Secondly, considering running parallel batch jobs across different R scripts using HPC then this approach of parallelsim provides a way for the optimal execution of those R scripts that use the available computer resources (e.g., computing capabilities) as efficient as possible. See the following bibliography:

- Eugster, M. J., Knaus, J., Porzelius, C., Schmidberger, M., & Vicedo, E. (2011). Hands-on tutorial for parallel computing with R. Computational Statistics, 26(2), 219-239.
- Knaus, J. Snowfall: Easier cluster computing (based on snow). R package version 1.84-6; 2013.
- [Parallel Computation in R.](https://bookdown.org/rdpeng/rprogdatascience/parallel-computation.html)

```R

## Examples 




```


On Extending R with C++ Capabilities: 
- Eddelbuettel, D., & Balamuta, J. J. (2018). [Extending R with C++: A brief introduction to Rcpp. The American Statistician, 72(1), 28-36.](https://www.tandfonline.com/doi/abs/10.1080/00031305.2017.1375990)

In practise the understanding of some important computing tools and techniques can assist in improving algorithmic efficiency and reduce the execution time. A classic book which covers related topics is:

- Tenenbaum, A. M., & Augenstein, M. J. (1986). Data structures using Pascal. Prentice-Hall, Inc.

Lastly, once a statistical/econometrics project containing R code which can be compiled has matured to the point of sharing it with collaborators or using it within a parallel computing environment, then the ideal way forward is to embed the code within an R package (see,  [Eddelbuettel and Balamuta (2018)](https://www.tandfonline.com/doi/abs/10.1080/00031305.2017.1375990)).     

## Random Generating Numbers 

> Question: What is randomness? Can we define what a 'random number' is?

When running a simulation study is important to be able to replicate the same results on different R sessions (on the same or different computers) using the same coding procedure. Thus, it is possible to reproduce the results by seeding the RNG, essentially choosing where to start in a long sequence of randomly generated numbers (say uniformly distributed within a prespecified interval).      

```R

## Random number generators

> set.seed( seed, kind = NULL, normal.kind = NULL )

## Example: 

> set.seed(1234)
> runif(5)
[1] 0.1137034 0.6222994 0.6092747 0.6233794 0.8609154
> runif(5)
[1] 0.640310605 0.009495756 0.232550506 0.666083758 0.514251141
> set.seed(1234)
> runif(5)
[1] 0.1137034 0.6222994 0.6092747 0.6233794 0.8609154

```

Furthermore, it is convenient to have the seed set randomly every time a program is run, so that the results are independent each time. To do this, one option is to use information from the computer itself, such as the current time or the process number assigned to the program that is running. 


# [A]. Optimization Techniques and Dimensionality Reduction

## [A1]. Portfolio Optimization Techniques 

Consider $r_N$ to be the N-dimensional vector of stock returns at time t. In other words, this is considered to be a strictly stationary time series of an (N x 1) observed vector of returns (e.g., log returns) with population moments given by 

$$ \boldsymbol{\mu} = \mathbb{E}[ \boldsymbol{r}_t ], \ \ \text{where} \ \boldsymbol{\mu} \in \mathbb{R}^{N \times 1}.$$

representing the vector of expected returns for a given time period, e.g., a historical period or a full sample period $t=1,...,T$. 

$$\boldsymbol{\Sigma}= \mathsf{Cov}[ \boldsymbol{r}_t ] = \mathbb{E}[ \boldsymbol{r}_t \boldsymbol{r}_t^{\top} ] - \boldsymbol{\mu} \boldsymbol{\mu}^{\top}.$$
representing the (unconditional) covariance matrix (often called volatility matrix) of stock returns for the time period $t=1,...,T$. 

Moreover, let $R_{t+1}$ denote the continuous random variable which represents the portfolio returns at time t+1 which is given by
$$R_{t+1} = \sum_{i=1}^N w_{i,t} r_{i,t+1} \equiv \boldsymbol{w}^{\top} \boldsymbol{r}_{t+1}.$$

Then, the traditional minimum-variance portfolio allocation problem

$$\underset{ \boldsymbol{w} \in \mathbb{R}^N }{ \mathsf{arg min} } \ \boldsymbol{w}^{\top} \boldsymbol{\Sigma} \boldsymbol{w} \ \ \ \text{subject to} \ \ \ \boldsymbol{w}^{\top} \boldsymbol{1} = 1.$$

where $w = ( w_1,..., w_N )$ represents the weights put on different assets and 1 is a N-dimensional vector with all entries being 1. Then, the optimal vector of weights which has a closed-form solution is given by the following expression: 

$$w^* = \frac{ \boldsymbol{\Sigma}^{-1} \boldsymbol{1} }{ \boldsymbol{1}^{\top} \boldsymbol{\Sigma}^{-1} \boldsymbol{1} }.$$

Therefore, for high-dimensional portfolio optimization problems, a challenging problem is the robust estimation of the inverse covariance matrix (precision matrix). In practise since the sample covariance matrix is used as a statistic of the corresponding population covariance matrix, then the resulting 'plug-in' portfolio vector has been widely adopted. However a related research question of interest is: How well does such a portfolio perform? A first step in this direction is to consider evaluating the portfolio performance - especially when comparing different estimates, using for instance stochastic dominance tests and/or via other portfolio performance measures. A second step in tackling this aspect is to investigate the asymptotic properties of the portfolio using different estimators for the covariance matrix by deriving appropriate error bounds in relation to different optimization constraints. 

Furthermore, the mean-variance portfolio optimization problem is defined as below

$$ \boldsymbol{w}^* =  \underset{ \boldsymbol{w} \in \mathbb{R}^N  }{ \mathsf{arg min} } \ \boldsymbol{w}^{\top} \boldsymbol{\Sigma} \boldsymbol{w}  \ \ \ \text{subject to} \ \ \ \boldsymbol{w}^{\top} \boldsymbol{\mu} = \boldsymbol{\mu}_0.$$

The closed-form solution of the above investement strategy is given by

$$\boldsymbol{w}^* = \frac{ \boldsymbol{ \mu}_0 }{ \boldsymbol{\mu}^{\top} \boldsymbol{\Sigma}^{-1} \boldsymbol{\mu} } \boldsymbol{\Sigma}^{-1} \boldsymbol{\mu}.$$

Similarly we can also consider Expected Utility Strategies, in which we incorporate the expected utility of an investor and thus the degree of risk aversion within the objective function as shown below 

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

In practise, the additional condition of having the weight vector being between zero and one is necessary when constructing the optimization problem to avoid infeasible solutions especially in the case where the optimal portfolio problem is constructed based on small number of observations. 

- To implement the aformentioned procedure in R, we can employ a specific type of optimization algorithms called ['Genetic Algorithm'](https://rpubs.com/Argaadya/550805) (GA) which employs an iterative procedure to converge to an approximate solution given a prespecified number of iterations.  

- The main idea of the the 'Genetic Algorithm' is that it considers a fitness function. The algorithm begins with a population of N "individuals" (nodes or stocks in the financial network) with DNA or weights $w_i$ for i = 1,...,N. The fitness of each individual is evaluated based on a fitness function in a manner bearing certain analogies to the exchange of DNA in biological organisms. Heuristically, the algorithm converges for a large number of replications which results to the optimal vector of weights $\hat{w}_i$ (see, Halbert White, (1989)).  

- Notice that for convex problems any local optimum is also a global optimum. However, in non-convex problems this is might not be the case. For example, when incorporating the eigenvalue centrality in the optimization problem further complexities, with respect to the convergence of a suitable algirithm, occur which we will need to check. The LASSO shrinkage for portfolio optimization can be also employed.  

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
           lower = rep(0, N), upper = rep(1, N), 
           maxiter = 1000, run = 200, optim = TRUE)
  
  optimal.weights <- weights(GA@solution)
  return( optimal.weights )
  
}# end of function

```

## Remarks:

- Further investment strategies can be considered, however the optimization problem has additional computational burden. Overall, we consider the optimization problem under investigation to be a learning problem with the optimal weights being random draws from the parameter space of portfolio weights such as  

$$Q =  \underset{ w_1,...,w_N \geq 0 }{  \mathsf{sup} }  \bigg[ \phi(w_1,...,w_N) | w_1,...,w_N  \bigg].$$

- Notice that for the implementation of the GA algorithm in R the maximum number of iterations (e.g., set to 1000) is considered to be a stopping rule for the convergence rate of the optimization problem. However, depending from the investment strategy and the corresponding constraints the 'optimal number' of iterations that ensure convergence to the approximate optimal vector can vary.  

- An alternative way of obtaining the solution of the optimal portfolio choice problem (minimum-variance investement strategy, i.e.,  Global Minimum Variance Portflio (GMVP)) is to use the least squares projection method (see, Maillet et al. (2015)). 

Define with $\boldsymbol{\Theta} = \boldsymbol{I} - \frac{ \boldsymbol{1} }{ \boldsymbol{1}^{\top} }$  and $\boldsymbol{Q}$ an $n \times (N-1)$ matrix which has as its columns the $(N-1)$ non-zero eigenvalues of the matrix  $\mathbf{\Theta}$ and it satisfies the properties $\mathbf{Q}^{\top}\mathbf{1} = 0 \ \text{and} \  \mathbf{Q}^{\top} \mathbf{Q} = \mathbf{I}_{N-1}$, then the solution of the GMVP using least squares regression is given as below

$$\hat{\mathbf{w}}^* = \frac{\mathbf{1}}{n} - \mathbf{Q} \hat{\gamma}^{*}.$$

which is equivalent to the following shrinkage norm: 

$$\hat{\gamma}^{*} = \underset{ \mathbf{\gamma \in \mathbb{R} } }{\text{argmin}} (  \frac{\mathbf{1}}{n} - \mathbf{Q} \hat{\gamma} )^{\top} \mathbf{ \hat{\Sigma} } (  \frac{\mathbf{1}}{n} - \mathbf{Q} \hat{\gamma} ).$$

A proof of the above equivalence with respect to the closed-form solution of the GMVP can be found in Maillet, Tokpavi and Vaucher (2015).

- Notice that eigenvalue analysis can be employed to the covariance matrix which allows to decompose the matrix in accordance with the rank statistics of the eigenvalues. Various methodologies exist in the literature to obtain these statistics (see, Yu (1991)). Moreover, when the covariance matrix is constructed based on a factor structure then further methodologies have been developed to obtain the common factors as well as to bound the estimation error when estimating eigenvalues (see, Fan et al. (2020)). 

## References

- Katsouris, C. (2021a). Optimal Portfolio Choice and Stock Centrality for Tail Risk Events. [arXiv preprint:2112.12031](https://arxiv.org/abs/2112.12031)
- Maillet, B., Tokpavi, S., & Vaucher, B. (2015). Global minimum variance portfolio optimisation under some model risk: A robust regression-based approach. European Journal of Operational Research, 244(1), 289-299.
- Brandt, M. W., & Santa‐Clara, P. (2006). Dynamic portfolio selection by augmenting the asset space. The journal of Finance, 61(5), 2187-2217.
- Fan, J., Guo, J., & Zheng, S. (2020). Estimating number of factors by adjusted eigenvalues thresholding. Journal of the American Statistical Association, 1-10.
- Fan, J., Liao, Y., & Shi, X. (2015). Risks of large portfolios. Journal of Econometrics, 186(2), 367-387.
- Fan, J., Zhang, J., & Yu, K. (2012). Vast portfolio selection with gross-exposure constraints. Journal of the American Statistical Association, 107(498), 592-606.
- Markowitz, H.M. 1952. Portfolio selection. Journal of Finance 7(1), March, 77–91.
- Markowitz, H.M. 1956. The optimization of a quadratic function subject to linear constraints. Naval Research Logistics Quarterly 3, 111–33.
- Markowitz, H.M. 1959. Portfolio Selection: Efficient Diversification of Investments. New Haven: Yale University Press. Reprinted, New York: John Wiley and Sons, 1970.
- Senneret, M., Malevergne, Y., Abry, P., Perrin, G., & Jaffres, L. (2016). Covariance versus precision matrix estimation for efficient asset allocation. IEEE Journal of Selected Topics in Signal Processing, 10(6), 982-993.
- White, H. (1989). Learning in artificial neural networks: A statistical perspective. Neural computation, 1(4), 425-464.
- Yu, K. B. (1991). Recursive updating the eigenvalue decomposition of a covariance matrix. IEEE Transactions on Signal Processing, 39(5), 1136-1145.

## Further Reading

- Antell, J., & Vaihekoski, M. (2019). Expected and realized returns in conditional asset pricing models: A new testing approach. Journal of Empirical Finance, 52, 220-236.
- Avella-Medina, M., Battey, H. S., Fan, J., & Li, Q. (2018). Robust estimation of high-dimensional covariance and precision matrices. Biometrika, 105(2), 271-284.
- Chen, J., Dai, G., & Zhang, N. (2020). An application of sparse-group lasso regularization to equity portfolio optimization and sector selection. Annals of Operations Research, 284(1), 243-262.
- Chi, Y., Lu, Y. M., & Chen, Y. (2019). Nonconvex optimization meets low-rank matrix factorization: An overview. IEEE Transactions on Signal Processing, 67(20), 5239-5269.
- David, R., & David, S. M. (2015). Statistics and Data Analysis for Financial Engineering: with R examples.
- Ledoit, O., & Wolf, M. (2022). Quadratic shrinkage for large covariance matrices. Bernoulli, 28(3), 1519-1547.
- Tyler, D. E., & Yi, M. (2020). Lassoing eigenvalues. Biometrika, 107(2), 397-414.
- Zhao, Z., Ledoit, O., & Jiang, H. (2020). Risk reduction and efficiency increase in large portfolios: leverage and shrinkage. University of Zurich, Department of Economics, Working Paper, (328).


## [A2]. Testing for Stochastic Dominance

A suitable statistical methodology for testing the statistical superiority of one distribution over another is to test for first-order stochastic dominance between two distribution functions. Tests for Stochastic Dominance are usually nonparametric tests with many applications when comparing observationally equivalent structures. Under certain regularity conditions stochastic dominance tests can be employed to evaluate the stochastic dominance of portfolios and portfolio returns. In particular, Post (2003) provide an extensive examination of various stochastic dominance test within the portfolio allocation framework. Morover, the framework of Linton et al. (2014) provide a testing framework of stochastic dominance efficiency of a given portfolio and can be particular useful in evaluating alternative options of efficient frontiers. Furthermore, Second-Order Stochastic Dominance can be evaluated within the framework of expected utility maximization. 

The main idea of such a statistical testing procedure is the interest in examining stochastic monotonic relationships for economic and financial relationships. Consider for example, two continuous random variables X and Y both supported on the space [0,1]. Firstly, the following assumptions should hold for the statistical validity when testing for stochastic dominance:

(i).   Investors use an expected utility maximization problem along with the portfolio returns. 
(ii).  The asset returns are assumed to be serially independent and identically distributed with 

$$E[ \mathbf{x} ] = \mathbf{\mu} \ \ \ \text{and} \ \ E[ (\mathbf{x} - \mathbf{\mu}).(\mathbf{x} - \mathbf{\mu})^{\top}] = \mathbf{\Sigma}.$$ 

(iii). Investors are allowed to diversify between assets with a corresponding N-dimensional vector of portfolio weights.

Secondly, a testing hypothesis of a stochastic dominance between two alternative distribution of returns has the form given below. Under the null hypothesis it holds that 

$$H_0: F_{Y|X}(y | x_1 ) \geq F_{Y|X}(y | x_2) \ \forall \ y, x_1, x_2 \in (0,1) \ \text{with} \ x_1 \leq x_2.$$  

While under the alternative hypothesis the space of interest is given by

$$H_1 : F_{Y|X}(y | x_1 ) \leq F_{Y|X}(y | x_2) \ \text{for some} \ y, x_1, x_2 \in (0,1) \ \text{with} \ x_1 \leq x_2.$$


## Example 2

Therefore, in order to be able to evaluate the effectiveness of the portfolio optimization problem and the corresponding portfolio returns induced by a risk matrix such as the financial connectedness matrix proposed by [Katsouris (2021)](https://arxiv.org/abs/2112.12031), we might also be interested to examine the robustness of the test under the distribution invariance assumption, or under changes in the network topology of nodes in the graph.

<p align="center">

<img src="https://github.com/christiskatsouris/Computational-Econometrics-R/blob/main/Notes/graph.jpg" width="650"/>

</p>  

>  Density Plots of Portfolio Risk based on a forecasting scheme with a rolling window of fixed length, $w = 250$, and an out-of-sample sequence of length $n = 100$. The particular choice was made heuristically by considering that the aim of our study is for the rolling window to be able to capture periods of increased interdependencies in the financial newtwork. Related aspect to examine is the optimal choice of the bandwidth parameter but we leave this aspect for future research. 

<p align="center">

<img src="https://github.com/christiskatsouris/Computational-Econometrics-R/blob/main/Notes/eigenvector_graph.jpg" width="550"/>

</p>  

## Remarks:

- The Density Plots displayed in figure 1 above provide a visual representation of the distribution of portfolio risk (for portfolio risk in this example we refer to the defintion of the quadratic form given by Fan et al. (2015)), when considering excluding central nodes based on different centrality measures, although such 'differences' have to be verified using a formal statistical methodology (e.g., such as a stochastic dominance test). In particular, due to the small sample size that we employed to obtain the out-of-sample forecast sequences, these results should be interpreted with this fact into consideration. In other words, in order to improve the approximation to the corresponding asymptotic distribution of portfolio risk (using the covariance-type matrix proposed by [Katsouris (2021)](https://arxiv.org/abs/2112.12031), which is based on heavy-tailed and dependent elements), we should consider a larger sample size, using for instance a suitable parallel programming technique.

- A second related aspect is the construction of the covariance matrix, for instance under the assumption of Gaussianity data versus based on heavy-tailed time series. In addition, when considering the risk matrix proposed by [Katsouris (2021)](https://arxiv.org/abs/2112.12031) has 'additional' dependence features (e.g., dependence across rows and colummns) that will need further examination. Furthermore, the estimation procedure has significant higher computational complexity, especially within a rolling window scheme, when the covariance matrix is replaced with the Gamma matrix proposed by [Katsouris (2021)](https://arxiv.org/abs/2112.12031). The reason is that the particular quantile connectedness matrix has a pairwise regression-based covariance-type structure. Therefore, some form of parallelism is required to obtain robust estimation results. 

- Moreover, the presence of finite-sample bimodality (as seen from the figures above), we conjecture that is a phenomenon occuring due to the presence of nonlinear cointegrating 'random objects' - that is, we introduce the idea of 'random regression objects' since the elements of the covariance-type matrix proposed by [Katsouris (2021)](https://arxiv.org/abs/2112.12031) are constructed from systems of pairwise (nodewise) quantile predictive regressions, and thus the risk matrix has entries with highly dependent and heavy-tailed data. In particular, the properties of covariance matrices with heavy-tailed data are discussed in the literature by [Professor Richard A. Davis](http://www.stat.columbia.edu/~rdavis/technical_reports.html) and his co-authors. Although we do not propose yet a study 'Towards a Unified Framework for Random Object Cointegrating Regressions' (a research agenda topic), we initially focus on exploring further the aforementioned properties both in the studies of [Katsouris (2021a)](https://arxiv.org/abs/2112.12031) and Katsouris (2021b). Furthermore, in order to formally investigate whether the presence of bimodality is due to heavy-tailedness or due to highly dependent data, we will need to first reconsider the notion of 'network dependence', and then via asymptotic theory analysis and related simulations, evaluate which effect stochastically dominates. Heuristically speaking, we can examine the finite-sample behaviour of such nonlinear functionals based on specific strong memory parameters, however the framework proposed by [Katsouris (2021a)](https://arxiv.org/abs/2112.12031) is much more challenging when considering the asymptotic theory analysis.   

- According to the recent working paper of Wang and Phillips (2022) titled: "A General Limit Theory for Nonlinear Functionals of Nonstationary Time Series", the following remark is given:

> Such bimodality is known to arise with self-normalized statistics and t ratios in other contexts, especially in the presence of heavy tailed data where individual large draws can dominate both the numerator and the denominator in the ratio (see Logan et al (1972)[16] and Fiorio et al (2010)[10]). The explanation of the phenomena in the present setting is unrelated to heavy tails but is instead related to strong dependence in the data (Remark 2.3, Wang and Phillips (2022)).

- The asymptotic theory for distributions of ratios of random variables is indeed an interesting aspect which has important implications for Survey Sampling theory. Furthermore, a remark on bimodality and weak instrumentation in structural equation modeling is given by [Phillips P.C.B (2006)](https://www.cambridge.org/core/journals/econometric-theory/article/remark-on-bimodality-and-weak-instrumentation-in-structural-equation-estimation/2C674680596B807FA356474AAC9C13AD).

- Lastly, we point-out that the framework proposed in the studies of Katsouris (2021a), discusses some important ideas with applications to other related fields as well. Specifically, the proposed algorithmic procedure of node centrality penalization (e.g., sequential node removal) which induces reduction of the tail dependency in the graph - without changing its main structure, is also related to the so-called persistent data structure as well as to persistent homology from the fields of computational topology and topological data analysis. These aspects indeed provide a motivation to examining more formally related estimation methodologies from the computational econometrics field for large data problems (see, Section [B] of this teaching page).    

<p align="center">

<img src="https://github.com/christiskatsouris/Computational-Econometrics-R/blob/main/Notes/results1.jpg" width="650"/>

</p>  

<p align="center">

<img src="https://github.com/christiskatsouris/Computational-Econometrics-R/blob/main/Notes/results2.jpg" width="650"/>

</p>  


## References
- Chetverikov, D., Wilhelm, D., & Kim, D. (2021). An adaptive test of stochastic monotonicity. Econometric Theory, 37(3), 495-536.
- Chen, L. Y., & Lee, S. (2018). Exact computation of GMM estimators for instrumental variable quantile regression models. Journal of Applied Econometrics, 33(4), 553-567.
- Chou, P. H., & Zhou, G. (2006). Using Bootstrap to Test Portfolio Efficiency. Annals of Economics & Finance, 7(2).
- Fan, J., Liao, Y., & Shi, X. (2015). Risks of large portfolios. Journal of Econometrics, 186(2), 367-387.
- Firpo, S., Galvao, A. F., Pinto, C., Poirier, A., & Sanroman, G. (2021). GMM quantile regression. Journal of Econometrics.
- Katsouris, C. (2021a). Optimal Portfolio Choice and Stock Centrality for Tail Risk Events. [arXiv preprint:2112.12031](https://arxiv.org/abs/2112.12031)
- Katsouris, C. (2021b). A Graph Topology Measure for a Time Series Regression-based Covariance Matrix with Tail Estimates. University of Southampton. Working paper. 
- Linton, O., Post, T., and Whang, Y.-J. (2014). Testing for the stochastic dominance efficiency of a given portfolio. The Econometrics Journal, 17(2):S59-S74.
- Post, T. (2003). Empirical tests for stochastic dominance efficiency. The Journal of Finance, 58(5):1905-1931.
- Phillips, P.C.B. (2006). A remark on bimodality and weak instrumentation in structural equation estimation. Econometric Theory, 22(5), 947-960.
- Wang and Phillips (2022). A General Limit Theory for Nonlinear Functionals of Nonstationary Time Series. Cowles Foundation Discussion Paper No. 2336. 


## Self-Assesment Questions:

1. What are the main econometric challenges when optimizing large portfolios? Explain the related issues with the well-known term 'curse of dimensionality'.

2. Can we construct an optimization function that relates centrality measures with risk measures? In particular, an estimation exercise which can shed light on computational aspects related to such a convex optimization approach is to implement an out-of-sample forecasting scheme (e.g., with a moving window of fixed length) based on the GMVP where penalization of central nodes is applied sequentially by removing one-central node at a time from the graph. The above computational problem can be also done by varying the centrality measures (see, this [page](https://github.com/christiskatsouris/Network-Analysis-R) for definitions) as a robustness check.    
3. How does the dependence structure of the risk matrix affect the porfolio optimization problem. Think for example the variance-covariance matrix versus the correlation matrix. What about if we replace the variance-covariance matrix with a risk matrix that measures tail dependency and therefore it corresponds to a covariance-type matrix based on heavy-tailed data with possibe high dependence?  

## Assignment 1

Following the following steps construct a bootstrap resampling simulation study where the portfolio performance is evaluated using different constraints (thus changing the investment strategy - focusing on the GMVP and the min-var optimization problems).

- Step 1: Choose an appropriate optimization methodology and set-up the constraint optimization problem using necessary matrix algebra tools (e.g., for instance the pseudo-inverse method for inversting S can be employed). 
- Step 2: Simulate from the Multivariate normal distribution the elements of the covariance matrix.
- Step 3: Compare the convergence properties of the optimization algorithm and the performance of the optimal portfolio using (i) the variance-covariance matrix and (ii) using a GARCH-based covariance matrix estimation method. 
- Step 4: Repeat the procedure B times for the bootstrap step and obtain the empirical variance and MSE in order to compare the different optimization methods and the two covariance matrix estimation methodologies mentioned in Step 3.


## [A3]. Principal Component Analysis 

The increase in computational power and the capacity to compute the asset allocation for large portfolios has moved the literature towards high intensive statistical and computational methods which require different approaches in terms of calculating the optimal portfolio asset allocation problem. In particular a seminal econometic framework in this direction is proposed by Chamberlain, G. and Rothschild, M. (1983) who examine the properties of the first two moments of the joint distribution of stock returns as the number of assets goes to infinity. The proposed covariance matrix decomposition captures the notion that the stochastic structure of asset returns is determined by a small number of factors, while the remaining factors can be considered negligible.  Further important studies in the direction of econometric techniques for estimating moments for large datasets include the methodology of principal components analysis, which we study in this section, for estimating large covariance matrices (see, Anderson et al (1963)) as well as the factor models estimation of large dimensions (see, Bai (2003)). Overall, when considering high dimensional data environments, estimation and robust statistical inference can be affected by the curse of dimensionality. Thus, one of the first techniques proposed in the literature to deal with the aspects of dimensionality reduction is the principal component analysis (PCA).  

- Revision Material related to Factor Models can be found in the [Notes](https://github.com/christiskatsouris/Computational-Econometrics-R/tree/main/Notes) folder of this teaching page. 

```R

## Examples 


```


## References

- Anderson, T. W. (1963). Asymptotic theory for principal component analysis. The Annals of Mathematical Statistics, 34(1), 122-148.
- Bai, J. (2003). Inferential theory for factor models of large dimensions. Econometrica, 71(1), 135-171.
- Bai, J., & Li, K. (2012). Statistical analysis of factor models of high dimension. The Annals of Statistics, 40(1), 436-465.
- Chamberlain, G. (1983). Funds, factors, and diversification in arbitrage pricing models. Econometrica: Journal of the Econometric Society, 1305-1323.
- Chamberlain, G. and Rothschild, M. (1983). Arbitrage, factor structure, and mean-variance analysis on large asset markets. Econometrica, 51(5):1281–1304.
- Tipping, M. E., & Bishop, C. M. (1999). Probabilistic principal component analysis. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 61(3), 611-622.
- James, W., & Stein, C. (1961). Estimation with quadratic loss Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability, Volume 1: Contributions to the Theory of Statistics, Berkeley.
- Tikhonov, A. N. (1963). On the solution of ill-posed problems and the method of regularization. In Doklady Akademii Nauk (Vol. 151, No. 3, pp. 501-504). Russian Academy of Sciences.



# [B]. Iterative Simulation and Estimation Methodologies

There are various iterative simulation techniques depending on the statistical inference problem under examination. A revision on resampling and Monte Carlo Methods from the statistics literature can be found in the [Notes](https://github.com/christiskatsouris/Computational-Econometrics-R/tree/main/Notes) folder of this teaching page. In this Section, we focus on estimation and computational algorithms when considering advanced econometric models based on economic theory. Specifically, we discuss the Minimum Distance Estimation (MDE) method and the Generalized Method of Moments (GMM).


## [B1]. Minimum Distance Estimation Method

A particular class of econometric models has parameters which are time-varying. For example, consider the case when evaluating the portfolio performance within an out-of-sample framework in which case one has to estimate the covariance matrix in each moving window. Then, a suitable econometric model to capture these time-varying moment conditions is given by 
$$y_t = g_t ( b_0 ) + u_t, \ \ \ \ \ \ \ \ \  t = 1,2,...,n$$

where the response variable $y_t$ is an ( m x 1) vector of observable random functions of discrete time and gt is a vector of known functions which depend on a (p x 1) vector of unknown parameters whose true value is denoted by $b_0$. In general, gt if a function of a number of exogenous variables as well as b0 so that the function is time dependent. The last component of the model is the vector ut of additive stochastic disturbances. In this case, an appropriate estimation methodology is to employ the iterated minimum distance estimator (MDE). Furthermore, it has been proved by  Phillips (1976) that there is an equivalent relation between the MDE estimator and the quasi-maximum likelihood (QML) estimator when the model disturbances are serially independent.  

Therefore, we consider any vector g ( n ) in which minimizes the following quadratic form 

$$ R_n( b ) = n^{-1} \sum_{t=1}^n ( y_n - g_t (b) )^{\prime} S ( y_n - g_t (b) )$$

given the observations $( y_t : t = 1,...,n )$ and some positive definite matrix S, the above expression is called a MDE of $b_0$. Then, concentrating the likelihood function with respect to b, we find that bn (S) minimizes the following expression

$$\text{log det} \ \left( n^{-1} \sum_{t=1}^n \left( y_n - g_t (b) \right)^{\prime} \left( y_n - g_t (b) \right) \right).$$

## Example 3

Consider the MDE of an AR(1) model using R. 


```R

## Install required packages in R
install.packages("KoulMde")
library(KoulMde)

set.seed(1234)
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


## [B2]. Generalized Method of Moments (GMM)

An alternative estimation technique to MDE and MLE is the Method of Moments, where the main idea is the mapping of sample moments to population moments. Furthermore, in the case of Generalized Method of Moments (GMM) there are more moments than parameters and thus econometric identification can be achieved. In addition, GMM can be employed regardless of any distributional assumptions which makes it an attractive alternative to MLE especially when Gaussianity assumption is not required.     

Define with 

$$g_n(  X, \theta ) = \frac{1}{n}  \sum_{i=1}^n g( X_i; \theta ) \ \ \ \text{such that} \ \ \mathbb{E} [ g_n(  X, \theta_0 ) ] = 0.$$

Then, the corresponding GMM estimator is obtained based on the following optimization function 

$$\hat{\theta}_n = \underset{  \theta \in \Theta  }{ \mathsf{argmin}  }  \ \ [ g_n(  X, \theta ) ]^{\prime} W [ g_n(  X, \theta ) ] $$

In other words, the statistical geometry of the GMM estimator minimizes the distance of sample counterparts to popultation parameters based on the underline orthogonality conditions and the distance matrix W.  


## Example 4

```R

## Install required packages in R
install.packages("gmm")
library(gmm)

set.seed(1234)

## Example 1: Simple GMM iterative method for parameter estimation
g1 <- function(theta,x) 
{
  m1 <- (theta[1]-x)
  m2 <- (theta[2]^2 - (x - theta[1])^2)
  m3 <- x^3-theta[1]*(theta[1]^2+3*theta[2]^2)
  f <- cbind(m1,m2,m3)
  return(f) 
}

Dg <- function(tet,x)
{
  G <- matrix(c( 1,2*(-tet[1]+mean(x)), -3*tet[1]^2-3*tet[2]^2,0, 2*tet[2],
                 -6*tet[1]*tet[2]), nrow=3,ncol=2)
  return(G)
}

n  <- 200
x1 <- rnorm(n, mean = 4, sd = 2)

##################
## Output:
##################

> print(results <- gmm(g1,x1,c(mu = 0, sig = 0), grad = Dg))
Method
 twoStep 

Objective function value:  0.01107618 

    mu     sig  
4.1775  2.0662  

Convergence code =  0

> summary(results)

Call:
gmm(g = g1, x = x1, t0 = c(mu = 0, sig = 0), gradv = Dg)

Method:  twoStep 

Kernel:  Quadratic Spectral(with bw =  0.214 )

Coefficients:
     Estimate     Std. Error   t value      Pr(>|t|)   
mu    4.1775e+00   1.4023e-01   2.9791e+01  5.1026e-195
sig   2.0662e+00   9.4414e-02   2.1884e+01  3.6835e-106

J-Test: degrees of freedom is 1 
                J-test   P-value
Test E(g)=0:    2.21524  0.13665

Initial values of the coefficients
      mu      sig 
4.162791 2.085222 

#############
Information related to the numerical optimization
Convergence code =  0 
Function eval. =  61 
Gradian eval. =  NA 

```

A useful application of computational econometrics beyond its use in portfolio optimization and GMM estimation, is the use of these algorithms in estimating quantile regression models, which has been pioneered by the seminar work of Koenker Roger and Zhijie Xiao. Quantile regression models have plethora of applications, some of these include: for instance when the econometrician is interested in causal inference based on conditional quantile functions as well as when estimating risk measures. A revision note on risk measures can be found below.

<p align="center">

<img src="https://github.com/christiskatsouris/Computational-Econometrics-R/blob/main/Notes/risk_measures.jpg" width="725"/>

</p>  

<p align="center">

<img src="https://github.com/christiskatsouris/Computational-Econometrics-R/blob/main/Notes/risk_measures1.jpg" width="600"/>

</p>  

<p align="center">

<img src="https://github.com/christiskatsouris/Computational-Econometrics-R/blob/main/Notes/risk_measures2.jpg" width="600"/>

</p>  

> Notice that the problem of 'Statistical Elicitability' arises when we consider the point forecast of risk measures such as the Expected Shortfall, Value at Risk, and other related risk management measures. Some related discussion can be found in the literature. Expain what are the main challenges when estimating such conditional dependent risk measures as well as proposed solutions to ensure that econometric identification is optimally achieved? 

## Assignment 2

Consider the GMM estimation for the econometric specification that corresponds to the quantile regression (QR) estimator as given by the following optimization function: 
$$\hat{ \beta } ( \tau ) \in \underset{ \beta \in \mathbb{R}^{p} }{  \text{arg min} } \frac{1}{n} \sum_{t=1}^n \rho_{\tau} \left( Y_t - X_t^{\prime} b \right).$$
Noticee that equivalently, the QR estimator $\hat{ \beta } ( \tau )$ is also the generalized method of moments (GMM) estimator based on the unconditional moment restriction given by 
$$\mathbf{E} \bigg[ \bigg( \tau - \mathbf{1} \bigg( Y \leq X^{\prime} \hat{ \beta } ( \tau ) \bigg) \bigg) . X \bigg] = 0.$$

HINT: Notice that when the CQF is modelled via a linear to the regressors function, such that $Q( Y | X ) = X^{\prime} \beta ( \tau )$ or $F_Y ( X^{\prime} \beta ( \tau ) \big| X ) = \tau$, then the coefficient $\beta ( \tau )$ satisfies the conditional moment restriction given by
$$\mathbf{E} \bigg[ \tau - \mathbf{1} \bigg( Y \leq X^{\prime} \hat{ \beta } ( \tau ) \bigg) | X \bigg] = 0.$$

## References

- Andrews, D. W., & Lu, B. (2001). Consistent model and moment selection procedures for GMM estimation with application to dynamic panel data models. Journal of econometrics, 101(1), 123-164.
- Andrews, D. W. (1997). A stopping rule for the computation of generalized method of moments estimators. Econometrica: Journal of the Econometric Society, 913-931.
- Chesher, A., & Smith, R. J. (1997). Likelihood ratio specification tests. Econometrica: Journal of the Econometric Society, 627-646.
- Firpo, S., Galvao, A. F., Pinto, C., Poirier, A., & Sanroman, G. (2021). GMM quantile regression. Journal of Econometrics.
- Hansen, L. (1982): “Large Sample Properties of Generalized Method of Moments Estimators,” Econometrica.
- Hansen, B. E., & Lee, S. (2021). Inference for iterated GMM under misspecification. Econometrica, 89(3), 1419-1447.
- Lof, M. (2014). GMM Estimation with Non‐causal Instruments under Rational Expectations. Oxford Bulletin of Economics and Statistics, 76(2), 279-286.
- Park, S. K., Ahn, S. K., & Cho, S. (2011). Generalized method of moments estimation for cointegrated vector autoregressive models. Computational statistics & data analysis, 55(9), 2605-2618.
- Stock, J. H., & Wright, J. H. (2000). GMM with weak identification. Econometrica, 68(5), 1055-1096.


## [B3]. Expected Maximization Algorithm (EM)

## Example 5



```R
# Examples



```


## References

- Wu, C. J. (1983). On the convergence properties of the EM algorithm. The Annals of statistics, 95-103.


# Reading List

$\textbf{[1]}$ Abadir, K. M., & Magnus, J. R. (2005). Matrix algebra (Vol. 1). Cambridge University Press.

$\textbf{[2]}$ Hayashi, F. (2001). Econometrics. Princeton University Press.

$\textbf{[3]}$ White, H. (1996). Estimation, inference and specification analysis (No. 22). Cambridge University Press.

$\textbf{[4]}$ Belsley, D. A., & Kontoghiorghes, E. (Eds.). (2009). Handbook of Computational Econometrics. John Wiley & Sons.

$\textbf{[5]}$ Boyd, S., Boyd, S. P., & Vandenberghe, L. (2004). Convex optimization. Cambridge University Press.

## Advanced Topics:

$\textbf{[1]}$ Corbae, D., Stinchcombe, M. B., & Zeman, J. (2009). An introduction to Mathematical Analysis for Economic Theory and Econometrics. Princeton University Press.

$\textbf{[2]}$ Golub, G. H., & Van Loan, C. F. (2013). Matrix computations. JHU press.

$\textbf{[3]}$ Chen, X. (2007). Large Sample Sieve Estimation of Semi-nonparametric Models. Handbook of Econometrics, 6, 5549-5632.

$\textbf{[4]}$ Zaslavski, A. J. (2010). Optimization on Metric and Normed spaces (Vol. 44). Springer Science & Business Media.

$\textbf{[5]}$ Kushner, H., & Yin, G. G. (2003). Stochastic Approximation and Recursive Algorithms and Applications (Vol. 35). Springer Science & Business Media.


# Learning Outcomes

1. Understand the basic philosophy and tools of computational algorithms for solving convex problems (e.g., optimal portfolio allocation, asset pricing problems, GMM estimation).
2. Be able to apply computational algorithms to data-driven methodologies (e.g., use of the iterated GMM for inference purposes). 
3. Understand the basic properties of convergence rates for computational algorithms (e.g., for computational efficieny and asymptotic theory analysis purposes). 
4. Be able to apply statistical testing procedures (such as tests for stochastic dominance) to evaluate the performance of optimization schemes. 
5. Be able to implement estimation and optimization techniques to economic, finance and actuarial data.
6. Be able to use Statistical/Econometric Programming Software such as R, Matlab or Stata.  

# Disclaimer

The author (Christis G. Katsouris) declares no conflicts of interest.

The proposed Course Syllabus is currently under development and has not been officially undergone quality checks. All rights reserved. 

Any errors or omissions are the responsibility of the author. 

# Acknowledgments

The author has benefited by participating in workshops and training sessions related to High Performance Computing both at the University of Southampton as well as at University College London (UCL). 

Further mathematical/theoretical applications related to Numerical Analysis, Scientific Computing or High Performance Computing can be also found on the website of [Prof. Robert Scheichl](https://katana.iwr.uni-heidelberg.de/people/rob/).

# Historical Background

> Standing on the shoulders of giants.
> 
> $\textit{''If I have been able to see further, it was only because I stood on the shoulders of giants."}$
> $- \text{Isaac Newton, 1676}$ 

$\textbf{Jack Carl Kiefer}$ (25 January 1924 – 10 August 1981) was an American mathematical statistician at Cornell University (1952 to 1979) and the University of California Berkeley (1979 to 1981). His research interests included the optimal design of experiments, which was his major research area, as well as a wide variety of topics in mathematical statistics. Much of Kiefer's research was on the design of experiments; the American Statistician obituary calls him "undoubtedly the foremost worker in optimal experimental design". However, he also made significant contributions to other areas of statistics and optimization, including the introduction of golden section search (his master's thesis work) the Dvoretzky–Kiefer–Wolfowitz inequality and the Bahadur-Ghosh-Kiefer representation (with R. R. Bahadur and J. K. Ghosh) (Source: Wikipedia).

$\textbf{Raghu Raj Bahadur}$ (30 April 1924 – 7 June 1997) was an Indian statistician considered by peers to be "one of the architects of the modern theory of mathematical statistics". He received his doctorate from the University of North Carolina under Herbert Robbins in 1950 after which he joined University of Chicago. He worked as a research statistician at the Indian Statistical Institute in Calcutta from 1956 to 1961. Bahadur spent the remainder of his academic career in the University of Chicago. He published numerous papers and is best known for the concepts of "Bahadur efficiency" and the Bahadur–Ghosh–Kiefer representation (with J. K. Ghosh and Jack Kiefer). Bahadur also framed the Anderson–Bahadur algorithm[9] along with Theodore Wilbur Anderson which is used in statistics and engineering for solving binary classification problems when the underlying data have multivariate normal distributions with different covariance matrices (Source: Wikipedia).  

# How to Cite a Website

See: https://www.mendeley.com/guides/web-citation-guide/
