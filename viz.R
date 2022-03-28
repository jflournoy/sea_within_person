#' ---
#' title: "Within-person reliability"
#' output: html_document
#' ---


#'
#' >Nezlek, J. B. (2017). A practical guide to understanding reliability in
#' studies of within-person variability. Journal of Research in Personality, 69,
#' 149–155. https://doi.org/10.1016/j.jrp.2016.06.020
#' 

#'
#' \begin{aligned}
#' \text{Voxels: } y_{ijk} &= \pi_{0jk} + \epsilon_{ijk}\\
#' \text{Occasions: } \pi_{0jk} &= \beta_{00k} + \rho_{0jk}\\
#' \text{Persons: } \beta_{00k} &= \gamma_{000} + \nu_{00k}\\
#' \\
#' \text{Combined: } y_{ijk} &= \gamma_{000} + \nu_{00k} + \rho_{0jk} + \epsilon_{ijk}\\
#' \end{aligned}
#' 
#' Voxel-level reliability for a parcel:
#' 
#' $$
#' \alpha = \frac{\tau_\pi}{\tau_\pi + \frac{\sigma^2}{j}}
#' $$
#'

#'
#' where $\tau_\pi$ indicates the variance of the random component of the
#' parameter $\pi$, $\sigma^2$ is the variance of the error term
#' $\epsilon_{ijk}$, and $j$ is the number of voxels.
#'
#' Important caveat:
#'
#' >Perhaps more important, similar to the assumptions made when using
#' Cronbach’s alpha, the items are assumed to be parallel. According to Shrout
#' and Lane (2012) this means that “(a) Items must be replicate indicators of
#' the same underlying latent trait, (b) responses to the items must be
#' independent after adjusting for the level of the latent trait, and (c) the
#' magnitude of error associated with each item must be the same.”
#' 

library(data.table)
library(brms)
library(ggplot2)
library(argparse)
library(e1071)
library(patchwork)

# stat: skew 
skew <- function(x) {
  xdev <- x - mean(x)
  n <- length(x)
  r <- sum(xdev^3) / sum(xdev^2)^1.5
  return(r * sqrt(n) * (1 - 1/n)^1.5)
}

# color_scheme_set("blue")
# ppc_stat(y, yrep1, stat = "skew", binwidth = 0.01) + 
#   xlim(0, .6) + 
#   legend_none()

parser <- ArgumentParser(description='Viz within-person alpha')
parser$add_argument('--subsample', action = 'store_true', help = 'run the model with 15 randomly selected voxels (sampled within person-sessions)')
# args <- parser$parse_args(c('--moresamples'))
args <- parser$parse_args()

SUBSAMP <- args$subsample

if(SUBSAMP){
  suffix <- '-sub'  
} else {
  suffix <- ''
}

setDTthreads(5)

task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

if(is.na(task_id)){
  this_ROI <- 196 #TO BE SET BY SLURM ENV VAR
  message('Not running as slurm array, setting ROI to some value: ', this_ROI)
} else {
  this_ROI <- task_id
  message('Running as slurm array with task id=ROI: ', task_id, '=', this_ROI)
}

file_search_regex <- paste0('alpha-fearGTcalm-ROI', this_ROI, '(-more_\\d{2})*', suffix, '.rds')

roi_fits <- dir('fits', pattern = file_search_regex, full.names = TRUE)
fits <- lapply(roi_fits, readRDS)

cfit <- brms::combine_models(mlist = fits)

prior_summary(cfit)
summary(cfit)
#parnames(cfit)
brms::neff_ratio(cfit, pars = c('sd_id__Intercept', 
                                'sd_id:sess__Intercept', 
                                'b_Intercept'))
bayesplot::mcmc_neff(brms::neff_ratio(cfit, pars = c('sd_id__Intercept', 
                                                     'sd_id:sess__Intercept', 
                                                     'b_Intercept', 'sigma')))
bayesplot::mcmc_acf(cfit, lags = 50, pars = c('sd_id__Intercept', 
                                              'sd_id:sess__Intercept', 
                                              'b_Intercept', 'sigma'))
bayesplot::mcmc_rhat(brms::rhat(cfit))
bayesplot::mcmc_trace(cfit, pars = c('sd_id__Intercept', 
                                     'sd_id:sess__Intercept', 
                                     'b_Intercept'))
brms::pp_check(cfit, type = 'dens_overlay')
brms::pp_check(cfit, type = 'ecdf_overlay')
library(patchwork)
brms::pp_check(cfit, nsamples = 1000, type = 'stat', stat = 'skew') + 
brms::pp_check(cfit, nsamples = 1000, type = 'stat', stat = 'kurtosis')

np <- nuts_params(cfit)
sum(np$Value[np$Parameter == 'divergent__'])
fit_array <- as.array(cfit)
bayesplot::color_scheme_set("darkgray")
bayesplot::mcmc_parcoord(fit_array, np = np, pars = c('sd_id__Intercept', 
                                                      'sd_id:sess__Intercept', 
                                                      'b_Intercept'),
                         np_style = bayesplot::parcoord_style_np(div_color = "red", div_size = .5, div_alpha = 1))
bayesplot::mcmc_pairs(fit_array, np = np, pars = c('sd_id__Intercept', 
                                                   'sd_id:sess__Intercept', 
                                                   'b_Intercept'),
                      off_diag_fun = c("hex"),
                      # off_diag_args = list(size = 0.75),
                      condition = bayesplot::pairs_condition(),
                      # transform = list('sd_id__Intercept' = 'log',
                      #                  'sd_id:sess__Intercept' = 'log'),
                      np_style = bayesplot::pairs_style_np(div_size = 4))

length(np$Value[np$Parameter == 'divergent__'])
window <- range(which(np$Value[np$Parameter == 'divergent__'] == 1) %% 4000)

bayesplot::mcmc_trace(fit_array, pars = "sd_id__Intercept", np = np, window = window) +
  xlab("Post-warmup iteration")
bayesplot::mcmc_trace(fit_array, pars = "sd_id:sess__Intercept", np = np, window = window) +
  xlab("Post-warmup iteration")
bayesplot::mcmc_trace(fit_array, pars = "b_Intercept", np = np, window = window) +
  xlab("Post-warmup iteration")

lp <- log_posterior(cfit)
bayesplot::mcmc_nuts_divergence(np, lp)

bayesplot::mcmc_areas_ridges(fit_array, pars = c('sd_id__Intercept', 
                                                 'sd_id:sess__Intercept'), prob = .95, prob_outer = .99) +
bayesplot::mcmc_areas(fit_array, pars = c('sigma'), prob = .95, prob_outer = .99)

vc <- VarCorr(cfit, summary = FALSE)
# alpha_model_fit
tau_pi <- vc$`id:sess`$sd^2
tau_beta <- vc$id$sd^2
sigma2 <- vc$residual__$sd^2
nvoxel <- length(cfit$data$y)/dim(unique(cfit$data[, c('sess', 'id:sess')]))[[1]]
alpha_post <- tau_pi / (tau_pi + sigma2/nvoxel)
pvar <- tau_pi / (tau_pi + tau_beta + sigma2)
post_df <- data.frame(alpha = alpha_post[,1], prop_var = pvar[,1])

knitr::kable(posterior_summary(post_df, robust = TRUE), digits = 3)

bayesplot::mcmc_areas(post_df, pars = 'alpha', prob = .95, prob_outer = .99, point_est = 'median') + 
  bayesplot::mcmc_areas(post_df, pars = 'prop_var', prob = .95, prob_outer = .99, point_est = 'median')
  