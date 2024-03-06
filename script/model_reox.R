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
library(parallel)


#import data

hyp_reo <-
  list.files(
    path = "DATA/HYP_reox",
    pattern = "\\.xlsx",
    all.files = TRUE,
    full.names = TRUE
  )

listhyp <- lapply(hyp_reo, read_xlsx)
names(listhyp) <- tools::file_path_sans_ext(basename(hyp_reo))

td1hyp <- read_xlsx("DATA/td1_hyp.xlsx")
td1hyp <- as.data.frame(td1hyp)

bfr_reo <-
  list.files(
    path = "DATA/BFR_reox",
    pattern = "\\.xlsx",
    all.files = TRUE,
    full.names = TRUE
  )

listbf <- lapply(bfr_reo, read_xlsx)
names(listbf) <- tools::file_path_sans_ext(basename(bfr_reo))

td1bfr <- read_xlsx("DATA/td1_bfr.xlsx")
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

td1nor <- read_xlsx("DATA/td1_nor.xlsx")
td1nor <- as.data.frame(td1nor)

#data formatting

listhyp <- lapply(listhyp, setNames, c("time", "tsi"))

listhyp1 <- lapply(names(listhyp), function(names) {
  bob <- listhyp[[names]]

  # Find corresponding 'td1' value in 'td1bfr'
  td1_value <- td1hyp$td1[td1hyp$name == names] - 1
  filtdata <- filter(bob, time >= floor(td1_value))
  return(filtdata)
})

names(listhyp1) <- tools::file_path_sans_ext(basename(hyp_reo))

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

listnor1 <- lapply(names(listnor), function(names) {
  bob <- listnor[[names]]

  # Find corresponding 'td1' value in 'td1bfr'
  td1_value <- td1nor$td1[td1nor$name == names] - 1
  filtdata <- filter(bob, time >= floor(td1_value))
  return(filtdata)
})

names(listnor1) <- tools::file_path_sans_ext(basename(nor_reo))
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

#same hyp

output_list2 <- lapply(names(listhyp1), function(name) {
  bob <- listhyp1[[name]]

  # Find corresponding 'td1' value in 'td1bfr'
  extracted_value <- td1hyp$td1[td1hyp$name == name]

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

#same for nor

output_list3 <- lapply(names(listnor1), function(name) {
  bob <- listnor1[[name]]

  # Find corresponding 'td1' value in 'td1bfr'
  extracted_value <- td1nor$td1[td1nor$name == name]

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

td1 <- td1nor %>% slice(which(name == "S5_1NOR_raw")) %>% select(td1) %>% as.numeric()

tsi_base <- listnor1$S5_1NOR_raw %>% slice(which(time == floor(td1))) %>% select(tsi) %>% as.numeric()

st2 <- data.frame(
  tau1 = c(0, 60),
  tau2 =c(20, 200),
  A1 = c(0, 80),
  A2 = c(-50, 50),
  TD2 = c(90,160))

mod2 <- nls2(
  tsi ~ tsi_base + ifelse(time < td1, 0 ,A1 * (1 - exp(-(time - td1) / tau1))) + (ifelse(time < TD2, 0 , A2 * (1 - exp(
    -(time - TD2) / tau2
  )))),
  data = listnor1$S5_1NOR_raw,
  start = st2,
  algorithm = "random-search",
  control = nls.control(maxiter = 100000)
)

coef(mod2)

#plot

pred2 <- augment(mod2)

pred3 <- predict(mod2)

plotfit(mod2, smooth = T)

tss <- sum((listnor1$S5_1NOR_raw$tsi - mean(listnor1$S5_1NOR_raw$tsi))^2)

residual <- residuals(mod2)
rss <- sum(residual^2)

1-(rss/tss)

#calculate coef
TSS <- list()

TSS <- as.data.frame(list.rbind(lapply(listnor1,function(bib){sum((bib$tsi - mean(bib$tsi))^2)})))
colnames(TSS) <- "TSS"

RSS <- as.data.frame(list.rbind(lapply(output_list3, function(bob){
  residuals <- residuals(bob)
  RSS <- sum(residuals^2)
  })))

colnames(RSS) <- "RSS"

rsquared <- cbind(RSS,TSS)

rsquared <- mutate(rsquared, rsquared = 1-(RSS/TSS))

coeff <- as.data.frame(list.rbind(lapply(output_list3, function(bob){coef(bob)})))

finalbfr <- as.data.frame(cbind(coeff, rsquared))
finalhyp <- as.data.frame(cbind(coeff, rsquared))
finalnor <- as.data.frame(cbind(coeff, rsquared))

#calculata A'

nb_rowbfr <- as.data.frame(lapply(listbf1, function(bub){nrow(bub)}))
nb_rowbfr <- t(nb_rowbfr) %>% as.data.frame()
colnames(nb_rowbfr) <- "nbrow"
nb_rowbfr <- rownames_to_column(nb_rowbfr)
finalbfr <- rownames_to_column(finalbfr)
colnames(td1bfr)[5] <- "rowname"
finalbfr <- full_join(finalbfr, nb_rowbfr, by = c("rowname"))
finalbfr <- full_join(finalbfr, td1bfr, by = c("rowname"))
finalbfr <- finalbfr %>% mutate(A1_prime = A1*(1-exp((-nbrow - td1)/tau1)))
finalbfr <- finalbfr %>% mutate(A2_prime = A2*(1- exp((-nbrow - TD2)/tau2)))
finalbfr <- finalbfr %>% mutate(Atot = A1_prime+A2_prime)

#plot all model

preds1 <- lapply(output_list1, function(bob){augment(bob)})
names(preds1) <- tools::file_path_sans_ext(basename(bfr_reo))

preds2 <- lapply(output_list2, function(bob){augment(bob)})
names(preds2) <- tools::file_path_sans_ext(basename(hyp_reo))

preds3 <- lapply(output_list3, function(bob){augment(bob)})
names(preds3) <- tools::file_path_sans_ext(basename(nor_reo))

listplot1 <- mapply(
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
  bob = preds3,
  bib = listnor1, bab = finalnor$rsquared,
  SIMPLIFY = FALSE
)

bfrplot <- ggarrange(plotlist = listplot, labels = names(listbf1), font.label = list(size =8, color = "black"), vjust = 0.8)
hypplot <- ggarrange(plotlist = listplot, labels = names(listhyp1), font.label = list(size =8, color = "black"), vjust = 0.8)
norplot <- ggarrange(plotlist = listplot1, labels = names(listnor1), font.label = list(size =8, color = "black"), vjust = 0.8)

norplot
bfrplot
hypplot
