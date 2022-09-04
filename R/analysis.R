# install.packages("metafor")
# install.packages("clubSandwich")
# install.packages("robumeta")
library(tidyverse)
library(metafor)
library(clubSandwich)
library(esc)


# get data
data <- readxl::read_excel("data/dataset2.xlsx",skip = 1) %>% janitor::clean_names()

dat<- data %>% select(author,total_patient,
                total_implant,bone_parameter,
                location,outcome,r_value,mean_age) %>%
  filter(author!= "Farré-Pagès (2011)") %>%
  mutate(est_id= paste0(1:n(),"_id"),
         b_para= case_when(
           bone_parameter=="Bone density (HU)"~ 1,
           bone_parameter=="Bone density (GV)"~ 2,
           bone_parameter=="Cortical thickness (mm)"~ 3),
         b_para= factor(b_para,levels = 1:3,labels = c("HU","GV","CT")),
         imp_per_pt= total_implant/total_patient,
         year= parse_number(author)) %>%
  group_by(b_para) %>%
  arrange(year,.by_group = T) %>%
  escalc(measure = "ZCOR",
         data = .,
         ri = r_value,
         ni = total_implant,
         slab = author)


mlabfun <- function(text, res) {
  list(bquote(paste(.(text),
                    " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                    ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                    I^2, " = ", .(formatC(res$I2, digits=1, format="f")), "%, ",
                    tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")")))}

textfun <- function(text, res) {
  bquote(paste(.(text),
               " (Q = ", .(formatC(res$QE, digits=2, format="f")),
               ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
               sigma^2.1, " = ", .(formatC(res$sigma2[1], digits=2, format="f")), "; ",
               sigma^2.2, " = ", .(formatC(res$sigma2[2], digits=2, format="f")), ";",
               " p(est) ", .(metafor:::.pval(res$pval, digits=3, showeq=TRUE, sep=" ")), ")",
  ))}


get_intcpt <- function(x){

  res<- conf_int(x,vcov = "CR2")
  res[1,]
}

results<- dat %>%
  group_by(outcome) %>%
  nest() %>%
  mutate(cov_mat= map(.x=data, ~with(.x,
                                     impute_covariance_matrix(vi=vi,
                                                              cluster = author,
                                                              r= 0.8)))) %>%
  mutate(intercept= map(.x=data,
                              .y=cov_mat,
                              .f = ~rma.mv(yi = yi~1,
                                           V = .y,
                                           slab = author,
                                           data = .x,
                                           random = ~ 1 | author/est_id,
                                           sparse = TRUE))) %>%

  mutate(bone= map(.x=data,
                              .y=cov_mat,
                              .f = ~rma.mv(yi = yi~1+b_para,
                                           V = .y,
                                           slab = author,
                                           data = .x,
                                           random = ~ 1 | author/est_id,
                                           sparse = TRUE))) %>%

  mutate(imp_pt= map(.x=data,
                              .y=cov_mat,
                              .f = ~rma.mv(yi = yi~1+imp_per_pt,
                                           V = .y,
                                           slab = author,
                                           data = .x,
                                           random = ~ 1 | author/est_id,
                                           sparse = TRUE))) %>%

  mutate(imp_pt_bone= map(.x=data,
                              .y=cov_mat,
                              .f = ~rma.mv(yi = yi~1+imp_per_pt+b_para,
                                           V = .y,
                                           slab = author,
                                           data = .x,
                                           random = ~ 1 | author/est_id,
                                           sparse = TRUE))) %>%
  mutate(across(contains(c("bone","imp")),
                ~map(.x, ~get_intcpt (.x)),.names ="{.col}_int" ))


# intercept only results-----
isq_intrecept <- results$intercept[[1]]
ptv_intrecept <- results$intercept[[2]]
itv_intrecept <- results$intercept[[3]]

# intercept + bone parameter as a moderator
isq_int_bone <- results$bone_int[[1]]
ptv_int_bone <- results$bone_int[[2]]
itv_int_bone <- results$bone_int[[3]]

# intercept + implant per patient as a moderator
isq_int_imp <- results$imp_pt_int[[1]]
ptv_int_imp <- results$imp_pt_int[[2]]
itv_int_imp <- results$imp_pt_int[[3]]

# access p value= isq_intrecept$QEp



# intercept + implant per patient + bone parameter as a moderator
isq_int_imp_bone <- results$imp_pt_bone_int[[1]]
ptv_int_imp_bone <- results$imp_pt_bone_int[[2]]
itv_int_imp_bone <- results$imp_pt_bone_int[[3]]


# forest plot ISQ ----------

isq<- results$data[[1]]
isq_ct <- isq %>% filter(b_para=="CT")
isq_gv <- isq %>% filter(b_para=="GV")
isq_hu<- isq %>% filter(b_para=="HU")

res_isq_ct <- rma(yi, vi, subset=(b_para=="CT"), data=isq)
res_isq_gv <- rma(yi, vi, subset=(b_para=="GV"),     data=isq)
res_isq_hu <- rma(yi, vi, subset=(b_para=="HU"),  data=isq)

xmin=-3
xmax=3.1
ymin=-10
ymax= nrow(isq)+ 15

forest(isq_intrecept,
       xlim = c(xmin,xmax),
       ylim= c(ymin,ymax),
       header = c("First author (Year)", "Correlation \n Coefficient [95% CI]") ,
       xlab = "Correlation Coefficient",
       psize = 1.1,
       slab = NA,
       cex=.8,
       order = b_para,
       rows = c(rev(24:37), rev(12:19), rev(3:7) ),
       # transf = convert_z2r,
       atransf = convert_z2r,
       ilab = cbind(location),
       ilab.xpos=c(-1.6),
       ilab.pos=4,
       mlab = NA,
       col = "blue",
       cex.lab = 0.9)
text(xmin, c(38.5,20.5,8.5), c("Bone density (HU)", "Bone density (GV)","Cortical thickness (CT)"
                        ), pos=4,cex = 0.9,font = 4)
text(xmin, rev(24:37), isq_hu$author  , pos=4, cex = 0.8)
text(xmin, rev(12:19), isq_gv$author  , pos=4, cex = 0.8)
text(xmin, rev(3:7), isq_ct$author  , pos=4, cex = 0.8)

text(-1.6, 38.5, "Location of the measure", cex = 0.75, font = 4,pos = 4 )
text(0,42, "Correlation between bone parameters and PIS (ISQ)", cex= 1)

addpoly(res_isq_hu, row=22.5,
        mlab=mlabfun("RE model for HU", res_isq_hu),
        col = "grey90",lty = "dashed")
addpoly(res_isq_gv, row= 10.5,
        mlab=mlabfun("RE model for GV", res_isq_gv),
        col = "grey90",lty = "dashed")
addpoly(res_isq_ct, row= 1.5,
        # mlab = NA,
        mlab=mlabfun("RE model for CT", res_isq_ct),
        col = "grey90",lty = "dashed")

text(-3,-1, "Overall multivariate model with RVE (only intercept)", cex = 0.75,pos = 4)
text(-3,-2, textfun("",isq_intrecept),cex = 0.6,pos = 4)

addpoly(x = isq_int_bone$beta,
        sei = isq_int_bone$SE,
        rows = -3,col = "grey",
        mlab = NA)
text(-3,-3, "Multivariate model with RVE (bone parameter)", cex = 0.75,pos = 4)
text(-3,-4, textfun("",results$bone[[1]]),cex = 0.6,pos = 4)

addpoly(x = isq_int_imp$beta,
        sei = isq_int_imp$SE,
        rows = -5,col = "grey",
        mlab = NA)
text(-3,-5, "Multivariate model with RVE (implants per patient)", cex = 0.75,pos = 4)
text(-3,-6, textfun("",results$imp_pt[[1]]),cex = 0.6,pos = 4)


addpoly(x = isq_int_imp_bone$beta,
        sei = isq_int_imp_bone$SE,
        rows = -7,col = "grey",
        mlab = NA)
text(-3,-7, "Multivariate model with RVE (bone parameter + implants per patient)", cex = 0.75,pos = 4)
text(-3,-8, textfun("",results$imp_pt_bone[[1]]),cex = 0.6,pos = 4)














