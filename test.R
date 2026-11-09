# Bounded synthetic checks only. No participant inputs or production Bootstrap.
args <- commandArgs(TRUE)
root <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else getwd()
source(file.path(root, "config.R"))
for (name in c("common.R", "data.R", "association.R", "mediation.R", "joint.R", "Table1.R"))
  sys.source(file.path(root, "analysis", name), envir = .GlobalEnv)
pass <- function(name) cat("PASS", name, "\n")

d <- data.frame(eid=1:6,cnt_older=c(0,1,2,0,0,2),cnt_brother=c(2,2,2,1,1,1),cnt_sister=0,
  adopted=c(0,0,0,1,-3,0))
z <- birth_order(d)
stopifnot(identical(z$birth_relative[1:3],c(0,.5,1)),is.na(z$birth_relative[4]),z$birth_relative[5]==0,is.na(z$birth_relative[6]))
pass("birth-position coding and eligibility")

set.seed(728); n <- 1000L
frame <- data.frame(eid=seq_len(n),birth_relative=runif(n),age=runif(n,40,70),sex=factor(sample(0:1,n,TRUE)),
  ethnicity=factor(sample(c("White","Mixed"),n,TRUE)),center=factor(sample(1:3,n,TRUE)),
  assessment_year=factor(sample(2007:2009,n,TRUE)),S=sample(2:10,n,TRUE))
y <- rbinom(n,1,plogis(-1+.2*frame$birth_relative+.01*frame$age))
s <- association_sample(frame,y,"birth_relative");got <- logistic_fit(s,"birth_relative")
expected <- summary(glm(s$formula,data=s$data,family=binomial()))$coefficients["birth_relative",]
stopifnot(abs(got$beta-expected[1])<1e-14,abs(got$SE-expected[2])<1e-14,abs(got$P-expected[4])<1e-14)
stopifnot(all(association_sample(frame,NULL,"birth_relative",cap=8,full=s)$data$S<=8))
pass("Logistic coefficient SE P and inherited S8 sample")

master <- data.frame(indicator_id=paste0("M",1:121),IHD_q_BH_global605_exploratory=.2,T2D_q_BH_global605_exploratory=.2,
  paired_Delta_IE_per_1000=1,paired_q_BH_within_contrast121_exploratory=.01,IHD_PM_percent=4,T2D_PM_percent=4)
master$IHD_q_BH_global605_exploratory[1:64]<-.01;master$IHD_PM_percent[1:44]<-6
selected <- select_mediators(master)
stopifnot(sum(selected$pass_three)==64,sum(selected$selected)==44)
master$IHD_PM_percent[44]<- -10
stopifnot(sum(select_mediators(master,check_counts=FALSE)$selected)==43)
pass("121 -> 64 -> 44 and signed PM rule")

effects <- mediation_effects(.05,.06,.07)
stopifnot(abs(effects$TE-effects$DE-effects$IE)<1e-12,abs(effects$PM-50)<1e-12)
a<-data.frame(bootstrap_id=1:3,IE=c(1,NA,3));b<-data.frame(bootstrap_id=2:4,IE=2:4)
p<-paired_effects(a,b);stopifnot(nrow(p)==1,p$bootstrap_id==3,p$difference==0)
stopifnot(identical(stars(c(.0005,.005,.02,.05,NA)),c("***","**","*","","")))
pass("risk units pairing and saved-q star thresholds")

C<-data.frame(x=seq_len(4));original<-C;weights<-rep(1,4);domains<-list(first="m1",second="m2")
draw<-function(b,C,uniforms)data.frame(m1=b+uniforms[,1],m2=.4*b+uniforms[,2])
predict<-function(b,M,C)cbind(IHD=plogis(-3+.2*b+.1*M$m1+.15*M$m2),T2D=plogis(-2-.2*b+.1*M$m1-.2*M$m2))
binding<-list(panel_id="fixture",bootstrap_id=0,sample_hash="s",population_hash="p",model_spec_hash="m",input_hash="i",rng_spec_hash="r")
before<-.Random.seed
v<-joint_product_stream(C,weights,domains,draw,predict,function(...)TRUE,uniform_dim=2,n_mc=2,seed=53,binding=binding)
stopifnot(identical(before,.Random.seed),identical(C,original))
for(e in c("IHD","T2D")) {
  r<-v$components[v$components$endpoint==e&v$components$allocation_order=="forward",]
  val<-function(k)r$effect[r$component_type==k]
  stopifnot(abs(val("TE")-val("DE")-val("IE_joint"))<1e-12,abs(sum(val("IIE"))+val("DEP")-val("IE_joint"))<1e-12)
}
pass("formal product-marginal decomposition closure and paired RNG preservation")

desc<-data.frame(birth_group=factor(c("First","First","Last"),levels=c("First","Middle","Last")),smoking=c("Never",NA,"Current"))
out<-describe_groups(desc,data.frame(column="smoking",label="Smoking",type="categorical"))
overall<-out[out$group=="Overall",];stopifnot(sum(overall$n)==3,abs(sum(overall$percent)-100)<1e-12)
pass("Table1 complete categories and missing denominator")

files<-list.files(root,pattern="[.]R$",full.names=TRUE,recursive=TRUE)
for(f in files)invisible(parse(f))
stopifnot(!any(grepl("(202[0-9]{5}|_v[0-9]|vendor|legacy|repairs)",substring(files,nchar(root)+2),ignore.case=TRUE)))
pass("flat code tree no historical worker directories")
