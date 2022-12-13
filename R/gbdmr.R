# beta is the methylation matrix, the row is the CpG name, and column is the participant
# phenotype is the covariable matrix, which should make all values numeric
# rho is the threshold for the block finding

gbdmr <- function(beta, phenotype, rho)
{
        library(maxLik)
        library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        library(minfi)
        library(geograbi)
        library(limma)
        library(corrplot)
        ###################################################################################################
        #This file contains source functions of gbdmr                                                 #
        #Given one CpG site, the basic model is: Z ~ beta0 + beta1*phenotype + beta2 * confounder1 + .... #
        ###################################################################################################
        #input:
        # Z: DNAm levels: n by L matrix; n is sample size; L is block size
        # X: covariates: n by p matrix; n is sample size; p is number of covariates; the first column must be phenotype
        #output
        # an maxLik object consisting of maximum likelihood, estimate, etc...
        get.mle.reg<-function(Z,X){
                X=as.matrix(X)
                Z=as.matrix(Z)
                p=ncol(X)-1
                n=nrow(Z)
                L=ncol(Z)
                log.Z=log(Z)
                log.1mZ=log(1-Z)
                log.denom=log(1+rowSums(Z/(1-Z)))
                starter.beta=mean((colMeans(Z)*(1-colMeans(Z))/apply(Z,MARGIN=2,var)-1)*(1-colMeans(Z)))
                if (p>0){
                        starter.alpha=as.matrix(lm((log.Z-log.1mZ) ~ X[,2:(p+1)])[[1]])
                        temp=rep(NA,length(p))
                        for (i in 2:(p+1)){
                                temp[i-1]=mean(starter.alpha[i,])
                        }
                        starter.alpha=c(starter.alpha[1,1:L],temp)
                        starter=c(starter.alpha,starter.beta)
                }else{
                        starter=c(as.vector(lm((log.Z-log.1mZ) ~ 1)[[1]]),starter.beta)
                }

                llk <- function(theta){
                        coef.intercept=theta[1:L]
                        if (p>0){
                                coef.covariate=theta[(L+1):(L+p)]
                                coef.beta=rbind(coef.intercept,matrix(coef.covariate,nrow=p,ncol=L))
                                beta=theta[p+L+1]
                                alpha=exp(X%*%coef.beta)*beta
                        }else{
                                coef.beta=matrix(coef.intercept,nrow=1)
                                beta=theta[p+L+1]
                                alpha=matrix(exp(coef.beta)*beta,nrow=n,ncol=L,byrow=TRUE)
                        }
                        if (L!=1){
                                nom1=sum(lgamma(rowSums(alpha)+ beta))
                                nom2=sum((alpha-1)*log.Z)-sum((alpha+1)*log.1mZ)
                                denom1=sum(lgamma(alpha))+n*lgamma(beta)
                                denom2=sum((rowSums(alpha)+beta)*log.denom)
                                return (nom1+nom2-denom1-denom2)
                        }else{
                                return (sum(log(dbeta(Z,alpha,rep(beta,length(alpha))))))
                        }
                }
                A=matrix(rep(0,p+L+1),nrow=1)
                A[p+L+1]=1
                B=as.matrix(0)
                res <- maxLik(llk, start = starter,constraints=list(ineqA=A, ineqB=B),iterlim=10000)
                return (res)

        }


        #output:joint p value of profile likelihood ratio/Wald method for that block (consiting of L CpG sites)
        ht<-function(Z,X.no.int){ #no intercept
                L=ncol(as.matrix(Z))
                X=cbind(1,X.no.int)
                res.constraint=get.mle.reg(Z,X[,-2]) #the first variable is the interesting phenotype
                res.unconstraint=get.mle.reg(Z,X)
                p.lr=(1-pchisq(2*(res.unconstraint$maximum-res.constraint$maximum),1)) #
                p=ncol(as.matrix(X))
                index=L+1
                Fisher=-(res.unconstraint$hessian[index,index])
                p.Wald=1-pchisq(res.unconstraint$estimate[index]^2*Fisher,1)
                return (list(p.Wald=p.Wald, p.lr=p.lr))

        }


        #############################################
        ###functions
        #############################################
        #input: x, y are both n by p matrices
        #output: a vector of length p indicating the column-wise correlation between x and y.

        colCors = function(x, y) {
                sqr = function(x) x*x
                if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
                        stop("Please supply two matrices of equal size.")
                x   = sweep(x, 2, colMeans(x))
                y   = sweep(y, 2, colMeans(y))
                cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
                return(cor)
        }

        #given a vector y, and a randomly generated vector as the same length of y
        #return a vector such that cor(y,output)=pho
        complement <- function(y, rho, x) {
                if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
                y.perp <- residuals(lm(x ~ y))
                return (rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2))
        }


        #cor(x,rowMeans(y))
        #logit function
        logit <- function(x) return (log(x/(1-x)))


        #################prepare the data
        beta <- beta
        t.beta <- t(beta)
        #get the map information for CpG sites
        map.data.subset <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        colnames(map.data.subset)
        #delete X and Y since our data delete the chr X and Chr Y
        map.data.subset = map.data.subset[map.data.subset$chr!='chrY'&map.data.subset$chr!='chrX',]
        map.data.subset <-map.data.subset$chr[map.data.subset$Name %in% rownames(beta)]
        table(map.data.subset)

        Z1 <- t.beta
        dim(Z1)
        ################################
        # get the block information ###
        #############################
        ##get start.loc.B and end.loc.B
        ##get start.loc.B and end.loc.B
        start.loc.B=c()
        end.loc.B=c()
        count=0
        for (chr.num in names(table((map.data.subset)))){
                count=count+1
                tbl=table((map.data.subset))

                which(names(tbl)==chr.num)
                index=(cumsum(tbl)[count]-tbl[count]+1):cumsum(tbl)[count]


                record_cor=colCors(Z1[,index[-length(index)]],Z1[,index[-1]])


                consecutive.table=rle(as.integer(record_cor > rho))
                block.len=consecutive.table[[1]]
                block.val=consecutive.table[[2]]

                temp.loc=1:length(index)
                temp.block.start=rep(0,sum(block.val==1,na.rm=T))
                temp.block.end=rep(0,sum(block.val==1,na.rm=T))
                for (i in 1:sum(block.val==1,na.rm = T)){
                        temp.block.start[i]=(cumsum(block.len)[block.val==1]-block.len[block.val==1]+1)[i]
                        temp.block.end[i]=cumsum(block.len)[block.val==1][i]+1

                }

                temp.start.loc.B=1:length(index)
                temp.end.loc.B=1:length(index)
                if (sum(block.val==1,na.rm = T)>0){
                        for (i in 1:sum(block.val==1,na.rm = T)){
                                temp.start.loc.B[temp.block.start[i]:temp.block.end[i]]=temp.block.start[i]
                                temp.end.loc.B[temp.block.start[i]:temp.block.end[i]]=temp.block.end[i]
                        }
                        temp.start.loc.B=unique(temp.start.loc.B)
                        temp.end.loc.B=unique(temp.end.loc.B)
                }


                if (count==1){
                        start.loc.B=temp.start.loc.B
                        end.loc.B=temp.end.loc.B
                }else{
                        start.loc.B=c(start.loc.B,temp.start.loc.B+cumsum(tbl)[count-1])
                        end.loc.B=c(end.loc.B,temp.end.loc.B+cumsum(tbl)[count-1])
                }

        }

        table(end.loc.B-start.loc.B+1)
        block_information <- table(end.loc.B-start.loc.B+1)
        total_compare_LR <- sum(table(end.loc.B-start.loc.B+1))

        L.vec <- as.numeric(names(table(end.loc.B-start.loc.B+1)))
        length(L.vec)
        #########################################################
        ########################################################
        # conduct the gbdmr test###

        n=nrow(t.beta)

        #number.positive.cpg.wald <- rep(NA,length(L.vec))
        number.positive.cpg.LR <- rep(NA,length(L.vec))
        count.L = 0
        index1 <- NA
        positi.record.LR1 <- NA
        start.loc.B1 <- NA
        end.loc.B1 <- NA
        for (L in L.vec) {
                count.L = count.L +1

                positi.record.LR=rep(NA,length(which(end.loc.B-start.loc.B+1==L)))
                count=0
                for (i in which(end.loc.B-start.loc.B+1==L)){
                        count=count+1
                        index=start.loc.B[i]:end.loc.B[i]

                        X = as.data.frame(phenotype)
                        X = as.matrix(X,n,ncol(phenotype)) #three means we have three variables
                        res=ht(t.beta[,index],X)
                        positi.record.LR[count]=res$p.lr
                        index2 <- index
                        index1 <- c(index1, index2)
                        positi.record.LR1[index[1]] <- positi.record.LR[count]

                }
                number.positive.cpg.LR[count.L] <- L * length(positi.record.LR[positi.record.LR < 0.05/total_compare_LR])
        }

        # combine the position information and p-value
        positi.record.LR2 <- positi.record.LR1[!is.na(positi.record.LR1)]
        length(positi.record.LR2)
        result_data_LR <- cbind(start.loc.B,end.loc.B,positi.record.LR2)
        result_data_LR <- as.data.frame(result_data_LR)
        colnames(result_data_LR) <- c("start","end","p_value")
        total_compare_LR
        # find the final DMRs
        DMR_LR <- result_data_LR[result_data_LR$p_value < 0.05/total_compare_LR,]
        DMR_LR
        DMR_LR$n <- DMR_LR$end - DMR_LR$start + 1
        total_cpg <- sum(DMR_LR$n)
        total_cpg

        #
        posit <- NA
        for (i in 1:nrow(DMR_LR)) {
                posit1 <- DMR_LR$start[i]:DMR_LR$end[i]
                posit <- c(posit,posit1)

        }
        posit <- posit[-1]


        #locate the cpg sites

        sig_cpg <- beta[posit,]
        cpg_name_LR <- c(rownames(sig_cpg))
        cpg_name_LR
        return(list("gbdmr" = DMR_LR, "positive_cpg" = cpg_name_LR, "gbdmrall" = result_data_LR))

}


























