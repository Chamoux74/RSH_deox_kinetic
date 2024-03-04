library(AgroReg)
library(pracma)
library(nlsr)
library(nlstools)
library(broom)
library(nls.multstart)
library(ggplot2)
library(tidyverse)
library(purrr)
library(rTPC)
library(nls2)
library(readxl)
library(rlist)
library(ggpubr)


#import data

hyp_reo <-
  list.files(
    path = "DATA/HYP_reox",
    pattern = "\\.xlsx",
    all.files = TRUE,
    full.names = TRUE
  )

listex <- lapply(hyp_reo, read_xlsx)
names(listex) <- tools::file_path_sans_ext(basename(hyp_reo))

bfr_reo <-
  list.files(
    path = "DATA/BFR_reox",
    pattern = "\\.xlsx",
    all.files = TRUE,
    full.names = TRUE
  )

listbf <- lapply(bfr_reo, read_xlsx)
names(listbf) <- tools::file_path_sans_ext(basename(bfr_reo))

td1bfr <- read_xlsx("DATA/td1_BFR.xlsx")
td1bfr <- as.data.frame(td1bfr)

nor_reo <-
  list.files(
    path = "DATA/NOR_reox",
    pattern = "\\.xlsx",
    all.files = TRUE,
    full.names = TRUE
  )

listnor <- lapply(nor_reo, read_xlsx)
names(listnor) <- tools::file_path_sans_ext(basename(nor_reo))

#data formatting

listex <- lapply(listex, setNames, c("time", "tsi"))

listex <- lapply(listex, function(bob){filter(bob, time >= -3)})

listbf <- lapply(listbf, setNames, c("time", "tsi"))

listbf1 <- lapply(names(listbf), function(names) {
  bob <- listbf[[names]]

  # Find corresponding 'td1' value in 'td1bfr'
  td1_value <- td1bfr$td1[td1bfr$name == names] - 1
  filtdata <- filter(bob, time >= floor(td1_value))
  return(filtdata)
})

names(listbf1) <- tools::file_path_sans_ext(basename(bfr_reo))

listnor <- lapply(listnor, setNames, c("time", "tsi"))

listnor <- lapply(listnor, function(bob){filter(bob, time >= -3)})
##loop model

set.seed(123)

#biexp formula

biexp <- function(A1, tau1, tau2, A2, TD1, TD2, TSImin1) {
  biexp1 <-  TSImin1 + A1 * (1 - exp(-(time - TD1) / tau1))
  biexp2 <-  ifelse(time <= TD2, 0,A2 * (1 - exp(-(time - TD2) /tau2)))
  return(biexp1 + biexp2)
}

#dataframe containing parameters

st1 <- data.frame(
  tau1 = c(0, 60),
  tau2 =c(0, 200),
  A1 = c(0, 80),
  A2 = c(-50, 50),
  TD2 = c(80,250))

#lapply model

output_list1 <- list()


output_list1 <- lapply(names(listbf1), function(name) {
  bob <- listbf1[[name]]

  # Find corresponding 'td1' value in 'td1bfr'
  extracted_value <- td1bfr$td1[td1bfr$name == name]

  tsi_base <- bob %>% slice(which(time == floor(extracted_value))) %>% select(tsi) %>% as.numeric()

  #Fit the model using 'nls2'
  mod1 <- nls2(
    tsi ~ tsi_base + ifelse(time <= extracted_value, 0 ,A1 * (1 - exp(-(time - extracted_value) / tau1))) + (ifelse(time <= TD2, 0 , A2 * (1 - exp(
            -(time - TD2) / tau2
          )))),
    data = bob,
    start = st1,
    algorithm = "random-search",
    control = nls.control(maxiter = 100000)
  )

  return(mod1)
})

#model indivual for wrong fitting

td1 <- td1bfr %>% slice(which(name == "S8_1BFR_raw")) %>% select(td1) %>% as.numeric()

tsi_base <- listbf1$S8_1BFR_raw %>% slice(which(time == floor(td1))) %>% select(tsi) %>% as.numeric()

st2 <- data.frame(
  tau1 = c(0, 60),
  tau2 =c(5, 200),
  A1 = c(0, 80),
  A2 = c(-50, 50),
  TD2 = c(50,250))

mod2 <- nls2(
  tsi ~ tsi_base + ifelse(time < td1, 0 ,A1 * (1 - exp(-(time - td1) / tau1))) + (ifelse(time < TD2, 0 , A2 * (1 - exp(
    -(time - TD2) / tau2
  )))),
  data = listbf1$S8_1BFR_raw,
  start = st2,
  algorithm = "random-search",
  control = nls.control(maxiter = 100000)
)

coef(mod2)

#plot

pred2 <- augment(mod2)

pred3 <- predict(mod2)

plotfit(mod2, smooth = T)

ggplot(listbf1$S3_1BFR_raw) +
  geom_point(aes(time, tsi)) +
  geom_line(
    aes(pred2$time, pred2$.fitted),
    color = "red",
    linewidth = 1.5
  ) +
  geom_point(aes(x= -3, y= 61.56), color = "red", inherit.aes = T) +
  geom_point(aes(x= 2, y= 63.1), color = "purple", inherit.aes = T)

tss <- sum((listbf1$S3_1BFR_raw$tsi - mean(listbf1$S3_1BFR_raw$tsi))^2)

residual <- residuals(mod2)
rss <- sum(residual^2)

1-(rss/tss)

#calculate coef
TSS <- list()

TSS <- as.data.frame(list.rbind(lapply(listbf1,function(bib){sum((bib$tsi - mean(bib$tsi))^2)})))
colnames(TSS) <- "TSS"

RSS <- as.data.frame(list.rbind(lapply(output_list1, function(bob){
  residuals <- residuals(bob)
  RSS <- sum(residuals^2)
  })))

colnames(RSS) <- "RSS"

rsquared <- cbind(RSS,TSS)

rsquared <- mutate(rsquared, rsquared = 1-(RSS/TSS))

coeff <- as.data.frame(list.rbind(lapply(output_list1, function(bob){coef(bob)})))

finalbfr <- as.data.frame(cbind(coeff, rsquared))
finalhyp <- as.data.frame(cbind(coeff, rsquared))
finalnor <- as.data.frame(cbind(coeff, rsquared))

#plot all model

preds1 <- lapply(output_list1, function(bob){augment(bob)})
names(preds1) <- tools::file_path_sans_ext(basename(bfr_reo))

listplot <- mapply(
  function(bob, bib, bab) {
    plot <- ggplot(data = bib) +
      geom_point(aes(time, tsi)) +
      geom_line(
        aes(time, .fitted),
        data = bob,
        color = "red",
        linewidth = 1.5
      ) +
      annotate(
        "text",
        x = 180,
        y = -Inf,
        label = paste("R-squared:", round(bab,3)),
        color = "black", size = 3,
        vjust = -0.5
      )
    return(plot)
  },
  bob = preds1,
  bib = listbf1, bab = finalbfr$rsquared,
  SIMPLIFY = FALSE
)

bfrplot <- ggarrange(plotlist = listplot, labels = names(listbf1), font.label = list(size =8, color = "black"), vjust = 0.8)
hypplot <- ggarrange(plotlist = listplot, labels = names(listex), font.label = list(size =8, color = "black"), vjust = 0.8)
norplot <- ggarrange(plotlist = listplot, labels = names(listnor), font.label = list(size =8, color = "black"), vjust = 0.8)

norplot
bfrplot
hypplot
