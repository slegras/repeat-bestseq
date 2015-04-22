args <- commandArgs(trailingOnly=TRUE)

loadFile = function(path=path, file=file){
	data = read.table(paste(path, file, sep="/"), sep="\t", header=T);
	return(data)
}

loadFiles = function(path=path, prefix=prefix){
	files = list.files(path = path, pattern = prefix)
	data= data.frame()

	for (file in files){
		data = rbind(data, loadFile(path=path, file=file))
	}
	return(data)
}

runAnalysis = function(path=path, prefix=prefix){

	data.all = loadFiles(path=path, prefix=paste(prefix, "\\.[0-9]+\\.count", sep=""))

	data.all.sort = data.all[order(data.all$"Count", decreasing = T),]

	f = function(x){return(length(unlist(strsplit(toString(x), split=""))))}
	data.length = sapply(data.all.sort$"Sequence", f)

	data.full = data.frame(data.all.sort, length=data.length)

	png(paste(path,"/",prefix,"_MappingCount_vs_Length.png", sep=""), width = 600, height = 600)
	plot(data.full$"Count", data.full$"length", pch=20, cex=.5, xlab="Mapping count", ylab="Length of repeat")
	dev.off()

	write.table(data.full, file=paste(path,"/",prefix, "_final.tsv", sep=""), 
		sep="\t", col.names=T, row.names=F, quote=F)

}

runAnalysis(path=args[1], prefix=args[2])

#path="/ifs/illumina/slegras/S00000_Jachowicz/Answer1/tmp"
#prefix="bestSeq_2015-avril-21-183152"

