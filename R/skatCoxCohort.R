skatCoxCohort <- function(Z, formula, SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", data=parent.frame(), verbose =FALSE){
	require(survival)
	env <- environment()
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	}
    
    #fit null model: 
    nullmodel <- coxph(formula=formula,data=data)
    X<-model.matrix(nullmodel,data)
	rn<-row.names(model.frame(nullmodel,data=data))
	nullcoef<-coef(nullmodel)
    	
	##match snps in Z with master list in SNPInfo file 
	mysnps <- colnames(Z)
	
	SNPInfo[,aggregateBy] <- as.character(SNPInfo[,aggregateBy])
	which.snps.Z <- colnames(Z) %in% SNPInfo[,snpNames]
	which.snps <- match(mysnps[which.snps.Z],SNPInfo[,snpNames])
	
	nsnps <- sum(which.snps.Z)
	if(nsnps == 0){ 
		stop("no column names in Z match SNP names in the SNP Info file!")
	}
	n = nrow(Z)
	
	if(verbose){
    	cat("\n Calculating signed LRTs... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = nsnps, style = 3)
    	pb.i <- 0
    }
	
	##fit individual betas/se's
	maf <- numeric(nrow(SNPInfo))
	maf[na.omit(which.snps)] <- colMeans(Z,na.rm=TRUE)[which.snps.Z]/2
	
	zlrt <- numeric(nrow(SNPInfo))
	zlrt[na.omit(which.snps)] <- apply(Z[,which.snps.Z, drop = FALSE],2,function(z){
		if(any(is.na(z))){
			mz <- mean(z, na.rm=TRUE)
			z[is.na(z)] <- mz
		}
        if (verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(nsnps/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))
		}
		 model<-coxlr.fit(cbind(z,X), nullmodel$y, nullmodel$strata, NULL,
	       init=c(0,nullcoef),coxph.control(iter.max=100),NULL,"efron",rn)
		return(sign(coef(model)[1])*sqrt(2*diff(model$loglik)))
		})
	
	if(verbose) close(pb)
	
	
	#deal with monomorphic SNPs
	zlrt[maf == 0] <- 0
	
	#differentiate missing from monomorphic:
	maf[!(SNPInfo[,snpNames] %in% colnames(Z))] <- -1

	#split into genes
	zlrt 	<- 	split(zlrt, SNPInfo[,aggregateBy])
	maf 	<- 	split(maf, SNPInfo[,aggregateBy])
	
	ngenes <- length(unique(SNPInfo[,aggregateBy]))
	if(verbose){
    	cat("\n Calculating covariance... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = ngenes, style = 3)
    	pb.i <- 0
    }
    
	##get covariance matrices:
	re <- as.list(by(SNPInfo[,snpNames], SNPInfo[,aggregateBy],function(snp.names){
		inds <- match(snp.names,colnames(Z))
		mcov <- matrix(0,length(snp.names),length(snp.names))
		if(length(na.omit(inds)) > 0){
			Z0 <- as.matrix(Z[,na.omit(inds),drop=FALSE])
			if(any(is.na(Z0))) Z0 <- apply(Z0,2,function(z){
				mz <- mean(z, na.rm=TRUE)
				z[is.na(z)] <- mz
				z
			})
				mod1 <- coxlr.fit(cbind(Z0,X), nullmodel$y, nullmodel$strata, NULL,
	       init=c(rep(0,ncol(Z0)),nullcoef),coxph.control(iter.max=0),NULL,"efron",rn)
				mcov[!is.na(inds), !is.na(inds)] <- solve(mod1$var[1:ncol(Z0),1:ncol(Z0)])
		}
		rownames(mcov) <- colnames(mcov) <- snp.names
		if(verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(ngenes/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))		  
		}
		return(mcov)
	}),simplify = FALSE)
	
	##aggregate
	for(k in 1:length(re)){
		re[[k]] <- list("scores" = zlrt[[k]]*sqrt(diag(re[[k]])), "cov" = re[[k]], "n" =n, "maf" = maf[[k]], "sey" = 1) 
	}
	if(verbose) close(pb)
	
	class(re) <- "skatCohort"
	return(re)
}
