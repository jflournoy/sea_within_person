fear <- data.frame(MASS::mvrnorm(100, mu = rep(0,15), Sigma = diag(15)))
neutral <- data.frame(MASS::mvrnorm(100, mu = rep(0,15), Sigma = diag(15)))

#low alphas
psych::alpha(fear, warnings = FALSE)$total$raw_alpha
psych::alpha(neutral, warnings = FALSE)$total$raw_alpha

#add physiological noise
noise_Sigma <- matrix(rep(.6, 15*15), nrow = 15)
diag(noise_Sigma) <- 1
noise <- data.frame(MASS::mvrnorm(100, mu = rep(0,15), Sigma = noise_Sigma))
#high alpha
psych::alpha(noise, warnings = FALSE)$total$raw_alpha

fear_noise <- fear + noise
neutral_noise <- neutral + noise
psych::alpha(fear_noise, warnings = FALSE)$total$raw_alpha
psych::alpha(neutral_noise, warnings = FALSE)$total$raw_alpha

#contrast fear and neutral
fear_gt_neutral <- fear_noise - neutral_noise
#very low alpha
psych::alpha(fear_gt_neutral, warnings = FALSE)$total$raw_alpha

#faster noise
fast_noise_1 <- data.frame(MASS::mvrnorm(100, mu = rep(0,15), Sigma = noise_Sigma))
fast_noise_2 <- data.frame(MASS::mvrnorm(100, mu = rep(0,15), Sigma = noise_Sigma))
fear_noise <- fear + fast_noise_1 + noise 
neutral_noise <- neutral + fast_noise_2 + noise
psych::alpha(fear_noise, warnings = FALSE)$total$raw_alpha
psych::alpha(neutral_noise, warnings = FALSE)$total$raw_alpha

#contrast fear and neutral
fear_gt_neutral <- fear_noise - neutral_noise
#very low alpha
psych::alpha(fear_gt_neutral, warnings = FALSE)$total$raw_alpha
