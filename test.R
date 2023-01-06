y_rng <- powerHLA::get_multinomial_rng(
  K = 5, af = c(0.8, 0.1, 0.05, 0.025, 0.025), B = 500, N = 100)
y <- powerHLA::get_multinomial_rng(
  K = 5, af = c(0.8, 0.1, 0.05, 0.025, 0.025), B = 1, N = 10000)
y <- as.vector(y)
a <- rep(x = 0.1, times = 5)
names(a) <- paste0("K", 1:length(a))

o <- powerHLA::get_power_analysis(
  y = as.vector(y), y_rng = y_rng, a = a)

ggplot(data = o$d_diff)+
  geom_errorbar(aes(x = allele, y = mean, ymin = X2.5., ymax = X97.5.,
                 group = B), width = 0.05,
             position = position_dodge(width = 0.85), col = "gray")+
  geom_point(aes(x = allele, y = mean, group = B),
             position = position_dodge(width = 0.85))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "darkgray")+
  theme_bw()


ggplot(data = o$d_lor)+
  geom_errorbar(aes(x = allele, y = mean, ymin = X2.5., ymax = X97.5.,
                    group = B), width = 0.05,
                position = position_dodge(width = 0.85), col = "gray")+
  geom_point(aes(x = allele, y = mean, group = B),
             position = position_dodge(width = 0.85))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "darkgray")+
  theme_bw()



ggplot(data = o$d_or)+
  geom_errorbar(aes(x = allele, y = mean, ymin = X2.5., ymax = X97.5.,
                    group = B), width = 0.05,
                position = position_dodge(width = 0.85), col = "gray")+
  geom_point(aes(x = allele, y = mean, group = B),
             position = position_dodge(width = 0.85))+
  geom_hline(yintercept = 1, linetype = "dashed", col = "darkgray")+
  theme_bw()
