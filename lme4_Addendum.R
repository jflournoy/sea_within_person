library(lme4)
#alpha_model_form <- bf(y ~ 1 + (1 | id/sess))

fit <- lmer(y ~ 1 + (1 | id/sess), data = roi_d)
summary(fit)
vc <- VarCorr(fit)
tau_pi <- diag(vc$`sess:id`)[[1]]
sigma2 <- attr(vc, 'sc')^2
nvoxel <- max(roi_d$v)
sprintf('Alpha: %0.3f', alpha_est <- tau_pi / (tau_pi + sigma2/nvoxel))
sprintf('Prop var: %0.3f', pvar <- tau_pi / (sum(unlist(lapply(vc, diag))) + sigma2))
