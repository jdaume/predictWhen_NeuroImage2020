library(sjPlot)
library(ggeffects)
library(splines)
library(dplyr)
library(lme4)
library(ggplot2)

load("meg_roi.RData") # contains ITPC and power in top 20% channels
df <- meg_roi

# rereference and scale data ----
df$condition <- relevel(df$condition, ref = 3)
df$itpcz = scale(df$itpc) # normalize
df$powz = scale(df$pow) # normalize

# model data ----
m1 <- lmer(itpcz ~ condition * ns(time, 10) + (1 + time | subj), data = df)

m2 <- lmer(itpcz ~ condition * ns(time, 10) + powz + (1 + time | subj), data = df)

# p <- ggpredict(m1, c("time[all]","condition","pow.blcz"))
# plot(p)
p <- ggpredict(m2, c("time[all]","condition"))
center <-  attr(df$itpcz,"scaled:center")
scalez <-  attr(df$itpcz,"scaled:scale")
p$semh <- (p$std.error+p$predicted)*scalez+center
p$seml <- (p$predicted-p$std.error)*scalez+center
p$predicted <- p$predicted*scalez+center
p$conf.low <- p$conf.low*scalez+center
p$conf.high <- p$conf.high*scalez+center
plot(p , 
     colors = c( "#666666","#20B2AA", "#663399"),
     show.title = FALSE,
     show.legend = FALSE,
     alpha = .8,
     use.theme = FALSE,
     show.y.title = FALSE
     )
theme_set(theme_ggeffects(base_size = 18, base_family = "Arial"))

anova(m1,m2)

# kick out other library before doing this!
library(lmerTest) # provides p values 
summary(m2) #gives model output with estimated df and p values using Satterthwaite (same p values as Kenward-Roger but much faster)
tab_model(m1,m2,show.ci = FALSE, show.stat = TRUE, string.stat = "t")

# library(parameters)
# p_value(m1,method="kenward")

