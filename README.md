MeiHMM: **M**is-segregation **E**rror **I**dentification through **H**idden **M**arkov **M**odels

MeiHMM segments trisomic chromosome 21 into blocks of two or three haplotype blocks and infers the stage of chromosome 21 nondisjunction in gametogenesis. It also identifies meiotic crossover event on chromosome 21. MeiHMM uses genotype data of Down syndrome proband, without requiring data from the parents.

MeiHMM is an R package and was built and tested on R 4.4.0. 

*Dependencies*

R package HMM (>= 1.0.1)

**Installation**

In **R**, run

	devtools::install_github("jjyanglab/MeiHMM")

**How to run**

When running MeiHMM for the first time, you need to download the reference files.

	download.file("https://github.com/jjyanglab/MeiHMM/releases/download/0.1.0/chr21.haplotypes.1000genomes.rda", destfile = "chr21.haplotypes.1000genomes.rda")
	download.file("https://github.com/jjyanglab/MeiHMM/releases/download/0.1.0/combined.af.data.rda", destfile = "combined.af.data.rda")
	load("chr21.haplotypes.1000genomes.rda")
	load("combined.af.data.rda")

To run a toy example (which takes 3-5 minutes)

	download.file("https://github.com/jjyanglab/MeiHMM/releases/download/0.1.0/example.txt", destfile = "example.txt")
	require(HMM)
	res <- MeiHMM("example.txt")
	plotMeiHMM(res)


