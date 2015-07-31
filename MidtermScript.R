####################################################
#
# R script for Midterm exam. 
# by Dan Jin. 
# last update: Oct. 15, 2012
# 
# To run this script, please do the following step in R.
#   1. source("MidtermScript.r")
#   2. result.ls=list() # it's better to sign the result to a list, since it will return large data.
#   3. result.ls=control.function()
#
# Result of this script:
#   1. Return one list, result.ls, that includes all results required.
#   2. Create 5 plots in png format. (See below for more details)
#
#
# Summary of the content in result.ls:
#   1. $beta.MLE: beta.MLE (note: includes all beta.MLE for each test);
#   2. $pval: p value (note: includes all p value for each test);
#   3. $hits.ls: all genotype marker (marker number and their p value) significant than Bonferroni correction;
#   4. $hit.set.ls: distinct sets;
#   5. $sig.hit.ls: the most significant marker in each set (marker number and the associated p value);
#   6. $sig.hit.ls: MAF for the most significant genotype markers in each set;
#   7. $sig.hit.GWAS: the most significant genotype marker in GWAS overall (marker number and its p value);
#   8. $corr.ls: correlation of Xa between this marker and one of the markers to its left and right colomn. correlation of Xa between this marker and marker 409
#
#Structure of result.ls:
# result.ls
# |-$beta.MLE
# |-$pval
# |-$hits.ls
# | |-$position
# | |-$pval
# |
# |-$hit.set.ls
# |-$sig.hit.ls
# | |-$position
# | |-$pval
# | |-$MAF
# |
# |-$sig.hit.GWAS
# | |-$position
# | |-$pval
# |
# |-$corr.ls
#   |-position.left
#   |-corr.left
#   |-position.right
#   |-corr.right
#   |-corr.409
#
#
# This script will also create 5 plots in png format.
#   1. Q2-Histogram of the phenotypes.png
#   2. Q3-beta.a_hat.png
#   3. Q3-beta.d_hat.png
#   4. Q5-Manhattan Plot.png
#   5. Q6-Manhattan Plot with Bonferroni correction.png
#
####################################################

#:::::::::::::::::Import & convert data:::::::::::::::::

import.sample <- function() {
	phenodata.df = read.table("midterm_phenotypes_fall12.txt", header = FALSE, sep = "\n")
	genodata.df = read.table("midterm_genotypes_fall12.txt", header = FALSE, sep = "")

	size = dim(genodata.df)
	sample.size = size[1]/2 # =500
	num.sample = size[2] # = 1194

	sample.ls = list()
	sample.ls$y = as.matrix(phenodata.df)
	sample.ls$g = genodata.df

	# convert genotype to Xa and Xd
	Xa = matrix(0, sample.size, num.sample)
	for (col in 1:num.sample) {
		unique.allel = unique(genodata.df[, col]) # find the unique allele
		allele.freq = c(sum(genodata.df[, col] == unique.allel[1]), sum(genodata.df[, col] == unique.allel[2])) # calculate allele frequence

		for (row in 1:sample.size) {
			if (genodata.df[2 * row - 1, col] == (genodata.df[2 * row, col])) {
				if (genodata.df[2 * row - 1, col] == unique.allel[match(max(allele.freq), allele.freq)]) {
					Xa[row, col] = 1
				} else if (genodata.df[2 * row - 1, col] == unique.allel[match(min(allele.freq), allele.freq)]) {
					Xa[row, col] = -1
				}
			}

		}
	}
	sample.ls$Xa = Xa
	sample.ls$Xd = 1 - 2 * abs(Xa)
	return(sample.ls)
}


#:::::::::::::::::F-statistics & p-value:::::::::::::::::

Fstats_pvals.ls <- function(sample.ls.y, sample.ls.Xa, sample.ls.Xd, sample.size) {
	#Create a list
	Fstats_pvals.ls = list()

	#Calculate F statistic
	##Calculate y bar
ybar = sum(sample.ls.y)/sample.size

	ones = rep(1, sample.size)
	x = cbind(ones, sample.ls.Xa, sample.ls.Xd)

	##Calculate beta.MLE
	Fstats_pvals.ls$beta.MLE = solve((t(x) %*% x)) %*% t(x) %*% sample.ls.y

	##Calculate y hat
	yhat = x %*% Fstats_pvals.ls$beta.MLE

	## SSM, SSE, MSM, MSE
	SSM = sum((yhat - ybar)^2)
	SSE = sum((sample.ls.y - yhat)^2)
	df.M = 2
	df.E = sample.size - 3
	MSM = SSM/df.M
	MSE = SSE/df.E
	Fstats_pvals.ls$Fstats <- MSM/MSE


	#Calculate p value
	Fstats_pvals.ls$pval <- pf(Fstats_pvals.ls$Fstats, df1 = df.M, df2 = df.E, lower.tail = FALSE, log.p = FALSE)
	Fstats_pvals.ls$log10pval <- -log10(Fstats_pvals.ls$pval)

	#Return
	return(Fstats_pvals.ls)
}


#:::::::::::::::::Control Function:::::::::::::::::

control.function <- function() {

	result.ls = list() # store all the results together in one list! Please refer to the top of this file for a summary and the structure of result.ls.


	#--------------Import sample--------------
	sample.ls = list()
	sample.ls = import.sample()
	genodata.df = sample.ls$g
	size = dim(sample.ls$Xa)
	sample.size = size[1]
	num.sample = size[2]


	#--------------Calculate f-test, p-value, beta.MLE--------------
	
	F_p.ls = list()
	beta.MLE = matrix(0, 3, num.sample)
	for (i in 1:num.sample) {
		Fstats_pvals.ls <- Fstats_pvals.ls(sample.ls$y, sample.ls$Xa[, i], sample.ls$Xd[, i], sample.size)
		F_p.ls$Fstats[i] = Fstats_pvals.ls$Fstats
		F_p.ls$pval[i] = Fstats_pvals.ls$pval
		F_p.ls$log10pval[i] = Fstats_pvals.ls$log10pval
		beta.MLE[, i] = Fstats_pvals.ls$beta.MLE
	}
	F_p.ls$beta.MLE = beta.MLE

	result.ls$beta.MLE = beta.MLE
	result.ls$pval = F_p.ls$pval
	result.ls$s1="############################################"


	#--------------Plotting--------------
	
	#QUESTION #2
	png("Q2-Histogram of the phenotypes.png", width = 500, height = 500)
	hist(sample.ls$y, main = "Histogram of the phenotypes", ylab = "Frequency", xlab = "Scaled height")
	dev.off()

	#QUESTION #3: Plot beta.a_hat and beta.d_hat
	position = seq(1, num.sample)
	png("Q3-beta.a_hat.png", width = 500, height = 500)
	plot(position, F_p.ls$beta.MLE[2, ], main = "beta.a hat vs. position", ylab = "beta.a hat", xlab = "position in genotype markers")
	dev.off()

	png("Q3-beta.d_hat.png", width = 500, height = 500)
	plot(position, F_p.ls$beta.MLE[3, ], main = "beta.d hat vs. position", ylab = "beta.a hat", xlab = "position in genotype markers")
	dev.off()

	#QUESTION #5: Manhattan plot
	png("Q5-Manhattan Plot.png", width = 500, height = 500)
	position = seq(1, num.sample)
	plot(position, F_p.ls$log10pval, main = "Manhattan Plot", ylab = "-log(p)", xlab = "position in genotype markers")
	dev.off()


	#QUESTION #6: Bonferroni correction
	alpha = 0.05
	bonferroni = alpha/num.sample
	bonferroni.log10 = -log10(bonferroni)
	position = seq(1, num.sample)
	y.bonferroni = rep(bonferroni.log10, num.sample)
	png("Q6-Manhattan Plot with Bonferroni correction.png", width = 500, height = 500)
	plot(position, F_p.ls$log10pval, main = "Manhattan Plot with Bonferroni Correction", ylab = "-log(p)", xlab = "position in genotype markers")
	lines(position, y.bonferroni, type = "l", col = "red")
	dev.off()


	#-------------finding hits--------------
	
	#QUESTION #7: List of hits
	alpha = 0.05
	bonferroni = alpha/num.sample
	bonferroni.log10 = -log10(bonferroni)

	hits.ls = list()
	h = 1
	for (p in 1:num.sample) {
		if (F_p.ls$log10pval[p] > bonferroni.log10) {
			hits.ls$positon[h] = p
			hits.ls$pval[h] = F_p.ls$log10pval[p]
#			hits.ls$hit.geno[h] = c(genodata.df[p])
			h = h + 1
		}
	}

	result.ls$hits.ls = hits.ls
	result.ls$s2="############################################"

	if (FALSE) {
		# How these hits distribute on the plot
		y = rep(0, 1, num.sample)
		y[hits.ls$positon] = 1
		png("test plot.png", width = 500, height = 500)
		plot(position, y)
		dev.off()
	} # end of if(FALSE)


	# Cluster hits into sets:
	## measure distance between two adjacent hits.
## set a threshold for clustering by calculating the mean distance. (this method depends on the following assumptions: 1.there should not be too many hits AND 2.the distance between two hits from two sets is much larger than the distance from two hits within a set)

	hit.position = c(hits.ls$positon)
	hit.dist = hit.position[-1] - hit.position[-length(hit.position)] # calculate the distance between two adjacent hits
	threshold = mean(hit.dist) # calculate the threshold of distance for clustering

	hit.set.ls = list() # list to store different sets of hits
	k = 1 # index for hit.set.ls
	hit.set = c() # vector to store the position of hits in the same set
	h = 1 # index for hit within a set
	hit.set[1] = hit.position[h]
	h = h + 1
	while (h <= length(hit.position)) {
		while (((hit.position[h] - hit.position[h - 1])) < threshold & (h <= length(hit.position))) {
			hit.set = c(hit.set, hit.position[h])
			h = h + 1
		}
		hit.set.ls[[k]] = hit.set
		k = k + 1
		hit.set = c() # empty his.set
		hit.set[1] = hit.position[h] # sign the current hit.position as the first element in the new hit.set
		h = h + 1
	}

	result.ls$hit.set.ls = hit.set.ls
	result.ls$s3="############################################"

	# Find the most significant marker in each sets:
	# QUESTION #8:Calculate MAF for each of the most significant marker in each sets:
num.set = length(hit.set.ls) # the number of distinct sets
	sig.hit.position = c()
	sig.hit.pval = c()
	sig.hit.MAF = c()
	for (i in 1:num.set) {
		# find the position with the most significant(smallest) p value in a set
		sig.pval.index = match(min(F_p.ls$pval[hit.set.ls[[i]]]), F_p.ls$pval[hit.set.ls[[i]]]) # index in hit.set.ls[[i]]
		sig.pval.position = hit.set.ls[[i]][sig.pval.index] # real marker position in genome

		sig.hit.position = c(sig.hit.position, sig.pval.position) # store the position of most sig p value from each set
		sig.hit.pval = c(sig.hit.pval, F_p.ls$pval[sig.pval.position]) # store the most sig p value from each set

		unique.allel = unique(genodata.df[, sig.pval.position]) # find the unique allele
		allele.freq = c(sum(genodata.df[, sig.pval.position] == unique.allel[1]), sum(genodata.df[, sig.pval.position] == 
			unique.allel[2])) # calculate allele frequence
		sig.hit.MAF = c(sig.hit.MAF, min(allele.freq)/sum(allele.freq)) # calculate MAF for the most sig hit of each sets
	}
	sig.hit.ls = list()
	sig.hit.ls$position = sig.hit.position
	sig.hit.ls$pval = sig.hit.pval
	sig.hit.ls$MAF = sig.hit.MAF

	result.ls$sig.hit.ls = sig.hit.ls
	result.ls$s4="############################################"


	# QUESTION #9: find the most significant marker in whole GWAS
	sig.hit.GWAS = list()
	sig.hit.GWAS$position = match(min(F_p.ls$pval), F_p.ls$pval) # marker position
	sig.hit.GWAS$pval = F_p.ls$pval[sig.hit.GWAS$position] # p value

	result.ls$sig.hit.GWAS = sig.hit.GWAS
	result.ls$s5="############################################"


	# Calculate correlations:
	corr.ls = list()
	# Correlation between Xa and Xa(column to its left)
	corr.ls$position.left = sig.hit.GWAS$position - 70
	corr.ls$corr.left = cor(sample.ls$Xa[, sig.hit.GWAS$position], sample.ls$Xa[, sig.hit.GWAS$position - 70])

	# Correlation between Xa and Xa(column to its right)
	corr.ls$position.right = sig.hit.GWAS$position + 3
	corr.ls$corr.right = cor(sample.ls$Xa[, sig.hit.GWAS$position], sample.ls$Xa[, sig.hit.GWAS$position + 3])

	# Correlation between Xa and Xa(marker 409)
	corr.ls$corr.409 = cor(sample.ls$Xa[, sig.hit.GWAS$position], sample.ls$Xa[, 409])

	result.ls$corr.ls = corr.ls
	
	return(result.ls)

}


