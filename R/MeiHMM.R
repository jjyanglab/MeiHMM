#' MeiHMM function
#'
#' MeiHMM segment chromosome 21 into blocks of two or three haplotypes
#' @param snp.file Path to the input file
#' @param ancestry The ancestry groups, can be a list in c("eur", "amr", "eas", "sas", "afr"). Use "tot" to ignore ancestry.
#' @return A list containing: segmentation, the segmentation results; positions, position of markers, af.doubled.allele, population allele frequency of two-copy alleles; oe.score, O/E ratio of hypothetical haplotypes
#' @export
MeiHMM <- function(snp.file, ancestry = c("tot")) {

	snp.info <- read.table(snp.file, sep ="\t")	
	pos <- as.numeric(snp.info[,2])
	
	xmin <- min(pos)
	xmax <- max(pos)
	
	vaf <- as.numeric(snp.info[,6]) / (snp.info[,5] + snp.info[,6])
	
	## q arm only, and exclude vaf = 1
	snp.info <- cbind.data.frame(snp.info, vaf)[pos > 1.2965e7 & !is.na(vaf) &vaf < 0.95, ]
	
	colnames(snp.info) <- c("chr", "position", "ref", "alt", "ref_count", "alt_count", "vaf")
	
	
	case.af.data <- combined.af.data[as.character(snp.info[, "position"]), ]
	
	
	population.afs <- sapply(ancestry, function(x) {
			
		
		total.population <- apply(case.af.data[, paste(x, c("A", "C", "G", "T"), sep = ".")], 1, sum)
		count.in.pop <- sapply(1:nrow(case.af.data), function(i) {
			case.af.data[i, paste0(x,  ".", snp.info[i, "alt"])]
			#dbsnp.data.here[as.character(snp.info[i, "position"]), paste0("SAMN10492705.", snp.info[i, "alt"])]
	
		})
		count.in.pop[count.in.pop == 0] <- 1
		population.af <- count.in.pop / total.population
	
	})
	
	population.afs <- apply(population.afs, 1, max)
	
	doubled.allele.af <- rep(NA, nrow(snp.info))
	alternative.allele.dose <- rep(NA, nrow(snp.info))
	

	# with ALT:REF:REF configuration
	ARR.msk <- snp.info[, "vaf"] > (1/3 - 0.2) & snp.info[, "vaf"] < (1/3 + gapping)
	doubled.allele.af[ARR.msk] <- (1-population.afs)[ARR.msk]
	# with ALT:ALT:REF configuration
	alternative.allele.dose[ARR.msk] <- 1
	AAR.msk <- snp.info[, "vaf"] > (2/3 - gapping) & snp.info[, "vaf"] < (2/3 + 0.2)
	doubled.allele.af[AAR.msk] <- population.afs[AAR.msk]
	alternative.allele.dose[AAR.msk] <- 2
	# low coverage positions
	low.cov.msk <- snp.info[, "ref_count"] + snp.info[, "alt_count"] < 30
	doubled.allele.af[low.cov.msk] <- NA
	alternative.allele.dose[low.cov.msk] <- NA
	
	to.use.msk <- (ARR.msk | AAR.msk) & !is.na(ARR.msk) & !is.na(AAR.msk)
	snp.info.to.use <- snp.info[to.use.msk, ]
	doubled.allele <- doubled.allele.af[to.use.msk]
	alternative.allele.dose <- alternative.allele.dose[to.use.msk]
	## arranging the two haplotypes
	haplotypes <- cbind(as.numeric(ARR.msk), as.numeric(AAR.msk))[to.use.msk, ]
	
	rownames(haplotypes) <- as.character(snp.info.to.use[, "position"])
	
	
	
	## haplotype analysis
	## not all SNPs appeared in 1000Genomes
	
	nmut.chr21 <- apply(chr21.haplotypes.1000genomes[, -(1:4)], 1, sum)
	pmut.chr21 <- nmut.chr21/ncol(chr21.haplotypes.1000genomes[, -(1:4)])
	names(pmut.chr21) <- rownames(chr21.haplotypes.1000genomes)

	
	msk.with.population.haplotypes <- rownames(haplotypes) %in% rownames(chr21.haplotypes.1000genomes)
	
	haplotypes.to.use <- haplotypes[msk.with.population.haplotypes,]
	
	population.haplotypes.here <- chr21.haplotypes.1000genomes[rownames(haplotypes.to.use), -(1:4)]
	
	
	firstbase <- c(t(population.haplotypes.here[1:(nrow(population.haplotypes.here) - 2), ]))
	
	middlebase <- c(t(population.haplotypes.here[2:(nrow(population.haplotypes.here) - 1), ]))
	
	lastbase <- c(t(population.haplotypes.here[3:(nrow(population.haplotypes.here)), ]))
	
	all.haplotypes <- paste0(firstbase, middlebase, lastbase)
	
	num.haplotypes.1000genome <- ncol(chr21.haplotypes.1000genomes) -4
	
	haplotype.counts <- sapply(1:(nrow(haplotypes.to.use) - 2), function(jj) {
	
		hap.counts <- rep(0, 8)
		names(hap.counts) <- c("000", "001", "010", "011", "100", "101", "110", "111")
		
		tb.counts <- table(all.haplotypes[(((jj - 1) * num.haplotypes.1000genome) + 1) : (jj * num.haplotypes.1000genome)])
		hap.counts[names(tb.counts)] <- tb.counts
		hap.counts
	})
	
	haplotype.counts <- t(haplotype.counts)
	rownames(haplotype.counts) <- rownames(haplotypes.to.use)[2:(nrow(haplotypes.to.use) - 1)]
	observed.haplotype.prop <- apply(haplotype.counts,2,"/", apply(haplotype.counts,1, sum))
	
	
	local.haplotypes <- cbind(paste0(haplotypes.to.use[1:(nrow(haplotypes.to.use) - 2), 1], haplotypes.to.use[2:(nrow(haplotypes.to.use) - 1), 1], haplotypes.to.use[3:(nrow(haplotypes.to.use)), 1]),
				paste0(haplotypes.to.use[1:(nrow(haplotypes.to.use) - 2), 2], haplotypes.to.use[2:(nrow(haplotypes.to.use) - 1), 2], haplotypes.to.use[3:(nrow(haplotypes.to.use)), 2]))
	
	
	#haps.na.mask <- apply(apply(actual.haps, 2, function(h) {grepl("NA", h)}), 1, sum) > 0
	
	paf.firstbase.mt <- pmut.chr21[rownames(haplotypes.to.use)[3:(nrow(haplotypes.to.use))]]
	paf.middlebase.mt <- pmut.chr21[rownames(haplotypes.to.use)[2:(nrow(haplotypes.to.use) - 1)]]
	paf.lastbase.mt <- pmut.chr21[rownames(haplotypes.to.use)[1:(nrow(haplotypes.to.use) - 2)]]
	
	
	paf.firstbase.wt <- 1 - pmut.chr21[rownames(haplotypes.to.use)[3:(nrow(haplotypes.to.use))]]
	paf.middlebase.wt <- 1 - pmut.chr21[rownames(haplotypes.to.use)[2:(nrow(haplotypes.to.use) - 1)]]
	paf.lastbase.wt <- 1 - pmut.chr21[rownames(haplotypes.to.use)[1:(nrow(haplotypes.to.use) - 2)]]
	
	
	expected.haplotype.prop <- cbind( 
			`000`=paf.firstbase.wt * paf.middlebase.wt * paf.lastbase.wt,
			`001`=paf.firstbase.wt * paf.middlebase.wt * paf.lastbase.mt,
			`010`=paf.firstbase.wt * paf.middlebase.mt * paf.lastbase.wt,
			`011`=paf.firstbase.wt * paf.middlebase.mt * paf.lastbase.mt,
			`100`=paf.firstbase.mt * paf.middlebase.wt * paf.lastbase.wt,
			`101`=paf.firstbase.mt * paf.middlebase.wt * paf.lastbase.mt,
			`110`=paf.firstbase.mt * paf.middlebase.mt * paf.lastbase.wt,
			`111`=paf.firstbase.mt * paf.middlebase.mt * paf.lastbase.mt)
	
	#actual.haps <- actual.haps[!haps.na.mask,]
	#hap.prop <- hap.prop[!haps.na.mask,]
	#observed.haplotype.prop <- observed.haplotype.prop[!haps.na.mask,]
	
	oe.score1 <- rep(NA, nrow(local.haplotypes))
	oe.score2 <- rep(NA, nrow(local.haplotypes))
	
	for(h in colnames(expected.haplotype.prop)) {
		oe.score1[local.haplotypes[,1] == h] <- observed.haplotype.prop[local.haplotypes[,1] == h, h]/expected.haplotype.prop[local.haplotypes[,1] == h, h]
		oe.score2[local.haplotypes[,2] == h] <- observed.haplotype.prop[local.haplotypes[,2] == h, h]/expected.haplotype.prop[local.haplotypes[,2] == h, h]
	}
	
	oe.score <- apply(cbind(oe.score1, oe.score2), 1, min)
	
	
	#oe.score[oe.score == 0] <- 0.0001
	
	observations <- rep("Others", nrow(snp.info.to.use))
	names(observations) <- as.character(snp.info.to.use[, "position"])
	observations[doubled.allele < rare_threshold] <- "Type1SNP"
	observations[rownames(observed.haplotype.prop)[oe.score < 0.1]] <- "Type2SNP"
	
	### HMM
	states <- c("2H", "3H")
	symbolsHMM <- c("Others", "Type1SNP", "Type2SNP")
	startProbs <- c(0.25, 0.75)   ## overall probability of 2H and 3H. Roghly 1/4 vs 1/4. actual value is 0.23379, estimated from first round 
	transProbs <- matrix(c(1-5e-10, 5e-10, 5e-10, 1-5e-10) , nrow = 2)  ## transition is rare but possible
	emissionProbs <- rbind(c(0.985, 0.0075, 0.0075), c(0.98, 0.000025, 1-0.98-0.000025))
	
	hmm <- initHMM(states, symbolsHMM, startProbs = startProbs, transProbs = transProbs, emissionProbs = emissionProbs)
	

	
	segmentation = viterbi(hmm, observations)
	
	oe.score.allpositions <- rep(NA, nrow(snp.info.to.use))
	names(oe.score.allpositions) <- as.character(snp.info.to.use[, "position"])
	oe.score.allpositions[rownames(observed.haplotype.prop)] <- oe.score
	
	return(list(segmentation = segmentation, positions = snp.info.to.use[, "position"], af.doubled.allele = doubled.allele,  oe.score = oe.score.allpositions))	
}




#' plotMeiHMM function
#'
#' plotMeiHMM generate figures for MeiHMM results for visual inspection
#' @param MeiHMM.results Returned list of function MeiHMM
#' @export

plotMeiHMM <- function(MeiHMM.results) {

	xmin <- min(MeiHMM.results$positions)
	xmax <- max(MeiHMM.results$positions)
	par(mar = c(0,4.16,4,1), mfrow = c(2,1), xpd=FALSE)
	
	plot(MeiHMM.results$positions, MeiHMM.results$af.doubled.allele, 
			log="y", #xlab = "chr21 position (million)", 
			ylab = "AF of doubled allele", xlim = c(xmin, xmax), ylim = c(1e-5, 1), pch=16, col = "grey20", cex = 0.5, axes = FALSE)
	#axis(1, at = (c(12:46)) * 1e6, labels = 12:46, las = 2)
	axis(2, at = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1), labels = c("1e-5", "1e-4", "1e-3", "1e-2", "1e-1", 1), las = 2)
	box()
	for (k in 12:46) {
		lines(c(1e6* (k + 0.5), 1e6* (k + 0.5)), c(1e-10, 100), lty = 2, col = adjustcolor("grey60", 0.8), lwd = 0.5)
	}
	
	changing.positions <- MeiHMM.results$positions[MeiHMM.results$segmentation != c(MeiHMM.results$segmentation[-1], "SomethingElse") | 
								MeiHMM.results$segmentation != c("SomethingElse", MeiHMM.results$segmentation[-length(MeiHMM.results$segmentation)])]
	changing.states <- MeiHMM.results$segmentation[MeiHMM.results$segmentation != c(MeiHMM.results$segmentation[-1], "SomethingElse") | 
								MeiHMM.results$segmentation != c("SomethingElse", MeiHMM.results$segmentation[-length(MeiHMM.results$segmentation)])]
	par(xpd = TRUE)
	for (i in 1:(length(changing.positions)/2)) {
		rect(changing.positions[2*i -1], 2, changing.positions[2*i], 2*5, col =  c("deepskyblue", "gold")[1+(changing.states[2*i] == "2H")], border = NA)
	}
	
	text(MeiHMM.results$positions[1]-1e6, 4, "#haplotype by HMM", adj = 1)
	
	par(mar = c(4.16,4.16,1,1), xpd=FALSE)
	
	plot(MeiHMM.results$positions, MeiHMM.results$oe.score, 
			log="y", xlab = "chr21 position (million)", 
			ylab = "score of hypothetical haplotype", xlim = c(xmin, xmax), ylim = c(1e-4, 5), pch=16, col = "grey20", cex = 0.5, axes = FALSE)
	axis(1, at = (c(12:46)) * 1e6, labels = 12:46, las = 2)
	axis(2, at = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = c("0", "0.001", "0.01", "0.1", "1"), las = 2)
	lines(c(-1, 1e8), c(0.1, 0.1), lty = 2, col = 2)
	box()
	for (k in 12:46) {
		lines(c(1e6* (k + 0.5), 1e6* (k + 0.5)), c(1e-10, 100), lty = 2, col = adjustcolor("grey60", 0.8), lwd = 0.5)
	}

}


