rm(list = ls())
library(TwoSampleMR)
library(MRPRESSO)
library(ieugwasr)
setwd('/media/desk15/iyun1417/MR')  
##############################################################################
#第一步数据处理
# 读取下载的数据文件
gwas_data <- readr::read_tsv("GCST90029028_buildGRCh37.tsv") 
gwas_data$N_CASE

# 转换为TwoSampleMR格式（需确保列名匹配）
exp <- format_data(
  gwas_data,
  type = "exposure",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  samplesize_col = 'N_CASE',
  pval_col = "p_value"
)

# 筛选p值并聚类
exp_filtered <-subset(exp, pval.exposure<5e-08)

#使用clump方法筛选合适的工具变量
exp_clumped <- clump_data(exp_filtered, clump_r2 = 0.001, clump_kb = 10000)


#输出结果
write.csv(exp_clumped, file="exposure.pvalue.csv", row.names=F)

#########################################################################################
#第二部去除F值
Ffilter=10        #F值过滤条件
inputFile="exposure.pvalue.csv"      #输入文件

#读取输入文件
dat<-read.csv(inputFile, header=T, sep=",", check.names=F)

#计算F检验值
dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))     #计算R2
dat$F<-dat$R2*(dat$samplesize.exposure-2)/(1-dat$R2)       #计算F检验值

#根据F值>10进行过滤, 删除弱工具变量
outTab=dat[dat$F>Ffilter,]
write.csv(dat, "exposure.F.csv", row.names=F)


###############################################################################
#第三步孟德尔随机化分析

#注释暴漏和结局--反向跑需要注意改
exposureName="Neuroticism"                      
outcomeName="IBS"   

exposureFile="exposure.F.csv"     
outcomeFile="GCST90016564_buildGRCh37.tsv"  
     
#读取除去F值暴露数据
exposure_dat=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                samplesize_col = "samplesize.exposure",
                                clump=FALSE)


#结局数据读取
gwas_data <- readr::read_tsv("GCST90029028_buildGRCh37.tsv") 

outcomeData<- format_data(
  gwas_data,
  type = "outcome",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  samplesize_col = 'n',
  pval_col = "p_value"
)

#合并相同SNP
outcomeTab<-merge(exposure_dat, outcomeData, by.x="SNP", by.y="SNP")
write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv")

#重新读入暴漏和结局数据
outcome_dat<-read_outcome_data(snps=exposure_dat$SNP,
                               filename="outcome.csv", sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta.outcome",
                               se_col = "se.outcome",
                               effect_allele_col = "effect_allele.outcome",
                               other_allele_col = "other_allele.outcome",
                               pval_col = "pval.outcome",
                               eaf_col = "eaf.outcome")

#注释暴漏和结局信息
exposure_dat$exposure=exposureName
outcome_dat$outcome=outcomeName
dat<-harmonise_data(exposure_dat=exposure_dat,
                    outcome_dat=outcome_dat)

#输出相同的SNP
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)

#多重检测
presso=run_mr_presso(dat)
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file="table.MR-PRESSO.csv")

#摘出重要方法
mrResult=mr(dat)

mr_method_list()$obj
mrResult=mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))

mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#7.异质性检验----MR-Egger P值
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#水平多效性【MR-Egger截距】------P值
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

#绘制图像

pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()


res_single=mr_singlesnp(dat)      
pdf(file="pic.forest.pdf", width=7, height=6.5)
mr_forest_plot(res_single)
dev.off()


pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()


pdf(file="pic.leaveoneout.pdf", width=7, height=6.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()

