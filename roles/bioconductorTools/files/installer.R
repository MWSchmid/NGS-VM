packageToInstall <- commandArgs(TRUE)[1]
source('http://bioconductor.org/biocLite.R')
if (packageToInstall == "bioconductor") {
	biocLite(ask=FALSE)
	print("Added")
} else {
	if (packageToInstall %in% installed.packages()[,"Package"]) {
		print('Already installed')
	} else {
		biocLite(packageToInstall, ask=FALSE)
		print("Added")		
	}
}
