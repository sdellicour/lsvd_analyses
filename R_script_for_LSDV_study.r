library(ape)
library(diagram)
library(fields)
library(igraph)
library(lubridate)
library(maptools)
library(MetBrewer)
library(RColorBrewer)
library(raster)
library(rgeos)
library(seraphim)
library(vegan)

# 1. Preparing the different alignment files for the global data set
# 2. Testing for a signal of recombination within the global alignment
# 3. Generating an haplotype network for the global data set (Network)
# 4. Inferring a maximum likelihood phylogeny for the global data set
# 5. Generating an overall sampling map for all the samples
# 6. Inferring the recombination breakpoints with GARD (www.datamonkey.org/GARD)
# 7. Preparing the input files for the analyses to perform in SPADS
# 8. Testing for a signal of recombination within the two wild type clades
# 9. Investigating the isolation-by-distance pattern based on the IID2 metric
# 10. Investigating the temporal signal and preparing the BEAST input files
# 11. Extracting the spatio-temporal information embedded in posterior trees
# 12. Mapping the dispersal history of LSVD lineages of the wild type clades

writingFiles = FALSE
savingPlots = FALSE

different_countries = c("Albania","Bulgaria","Croatia","Greece","Macedonia","Serbia",
						"Israel",
						"Turkey",
						"Kazachstan","Russia",
						"Bangladesh","China","India","Hong Kong","Taiwan","Thailand","Vietnam",
						"Kenya","Namibia","Nigeria","South Africa")
different_colours = c("#4676BB","#4676BB","#4676BB","#4676BB","#4676BB","#4676BB",
					  "#DE4327",
					  "#D1E5F0",
					  "#BF812D","#BF812D",
					  "#666666","#666666","#666666","#666666","#666666","#666666","#666666",
					  "#FAA521","#FAA521","#FAA521","#FAA521")

# 1. Preparing the different alignment files for the global data set

tab = read.csv("LSVD_all_alignment1.csv", head=T, sep=";")
txt1 = scan(paste0("LSVD_all_alignment1.fas"), what="", sep="\n", quiet=T)
tab[,1] = gsub(", complete genome","",tab[,1]); IDs = c(); txt2 = c()
for (i in 1:length(txt1))
	{
		if (grepl(">",txt1[i]))
			{
				txt2 = c(txt2, gsub(", complete genome","",txt1[i]))
			}	else	{
				if (grepl(">",txt2[length(txt2)]))
					{
						txt2 = c(txt2, txt1[i])
					}	else		{
						txt2[length(txt2)] = paste0(txt2[length(txt2)], txt1[i])
					}
			}
	}
IDs = txt2[grepl(">",txt2)]; IDs = gsub(">","",IDs)
collectionDates = as.character(dmy(gsub("\\.","-",tab[,"sampling_date"])))
for (i in 1:length(collectionDates))
	{
		if (is.na(collectionDates[i]))
			{
				date = tab[i,"sampling_date"]
				collectionDates[i] = as.character(dmy(gsub("\\.","-",date)))
			}
		if (is.na(collectionDates[i]))
			{
				date = tab[i,"sampling_date"]
				if (grepl("May/",date)) date = paste0(gsub("May/","",date),"-05")
				if (length(unlist(strsplit(date,"\\/"))) == 2)
					{
						date = paste0("15/",date)
						date = as.character(dmy(gsub("\\.","-",date)))
					}
				collectionDates[i] = date
			}
	}
for (i in 1:length(IDs))
	{
		index = which(tab[,"seqID"]==IDs[i])
		if (length(index) == 0) print(IDs[i])
		country = substr(tab[index,"country"],1,2)
		ID = unlist(strsplit(IDs[i],"\\."))[1]
		ID = substr(ID,nchar(ID)-3,nchar(ID))
		IDs[i] = paste0(ID,"-",country)
			# Network only allows 6 characters for the sequence IDs
	}
if (length(unique(IDs)) != length(IDs)) cat("Warning: some shorter sequence IDs are identical")
seqs1 = txt2[!grepl(">",txt2)]; seqs2 = matrix(nrow=length(seqs1), ncol=nchar(seqs1)[1])
for (i in 1:dim(seqs2)[1])
	{
		seqs2[i,] = unlist(strsplit(seqs1[i],""))
	}	# n.b.: if this loop creates continuous error messages, kill the AppleShed process with "kill -9 [PID]"
txt3 = paste0(length(IDs)," ",dim(seqs2)[2])
for (i in 1:length(IDs))
	{
		txt3 = c(txt3, paste0(IDs[i],"   ",paste(seqs2[i,],collapse="")))
	}
if (writingFiles) write(txt3, "LSVD_all_alignment1.phy")
seqs3 = seqs2
for (i in 1:dim(seqs3)[1])
	{
		if (seqs3[i,1] == "-") seqs3[i,1] = "£"
		if (seqs3[i,dim(seqs3)[2]] == "-") seqs3[i,dim(seqs3)[2]] = "£"
		indices = which(seqs3[i,] == "-")
		for (j in 1:length(indices))
			{
				if (seqs3[i,indices[j]-1] == "£") seqs3[i,indices[j]] = "£"
				if (seqs3[i,indices[j]+1] == "£") seqs3[i,indices[j]] = "£"
			}
	}
columnsToKeep = c()
for (i in 1:dim(seqs3)[2])
	{
		if (sum(seqs3[,i]!="£") == dim(seqs3)[1])
			{
				columnsToKeep = c(columnsToKeep, i)
			}
	}
seqs4 = seqs3[,columnsToKeep]
txt3 = paste0(length(IDs)," ",dim(seqs4)[2])
for (i in 1:length(IDs))
	{
		txt3 = c(txt3, paste0(IDs[i],"   ",paste(seqs4[i,],collapse="")))
	}
if (writingFiles) write(txt3, "LSVD_all_alignment2.phy")
	# "EU_verified_alignment1.phy" after 5'/3' trimming
columnsToKeep = c()
for (i in 1:dim(seqs2)[2])
	{
		if (sum(seqs2[,i]%in%c("a","c","t","g")) == dim(seqs2)[1])
			{
				columnsToKeep = c(columnsToKeep, i)
			}
	}
seqs3 = seqs2[,columnsToKeep]
txt3 = paste0(length(IDs)," ",dim(seqs3)[2])
for (i in 1:length(IDs))
	{
		txt3 = c(txt3, paste0(IDs[i],"   ",paste(seqs3[i,],collapse="")))
	}
if (writingFiles) write(txt3, "LSVD_all_alignment3.phy")
	# "EU_verified_alignment1.phy" only with the "a/c/t/g" sites
columnsToKeep = c()
for (i in 1:dim(seqs3)[2])
	{
		if (length(unique(seqs3[,i])) != 1)
			{
				columnsToKeep = c(columnsToKeep, i)
			}
	}
seqs4 = seqs3[,columnsToKeep]
txt4 = paste0(length(IDs)," ",dim(seqs4)[2])
for (i in 1:length(IDs))
	{
		txt4 = c(txt4, paste0(IDs[i],"   ",paste(seqs4[i,],collapse="")))
	}
if (writingFiles) write(txt4, "LSVD_all_alignment4.phy")
	# "EU_verified_alignment3.phy" only with the polymorphic sites
metadata = matrix(nrow=length(IDs), ncol=5); metadata[,1] = IDs
for (i in 1:length(IDs))
	{
		metadata[i,2] = unlist(strsplit(IDs[i],"-"))[length(unlist(strsplit(IDs[i],"-")))]
		index = which(substr(different_countries,1,2)==metadata[i,2])
		metadata[i,3] = different_colours[index]
	}
decimalDates = rep(NA, length(collectionDates))
for (i in 1:length(decimalDates))
	{
		if (length(unlist(strsplit(collectionDates[i],"-"))) == 3)
			{
				decimalDates[i] = decimal_date(ymd(collectionDates[i]))
			}
		if (length(unlist(strsplit(collectionDates[i],"-"))) == 2)
			{
				decimalDates[i] = decimal_date(ymd(paste0(collectionDates[i],"-15")))
			}
		if (length(unlist(strsplit(collectionDates[i],"-"))) == 1)
			{
				decimalDates[i] = as.numeric(collectionDates[i])+0.5
			}
	}
names(decimalDates) = IDs; colourScale = colorRampPalette(brewer.pal(11,"BrBG"))(101)
colIndices = (((decimalDates-min(decimalDates,na.rm=T))/(max(decimalDates,na.rm=T)-min(decimalDates,na.rm=T)))*100)+1
if (savingPlots)
	{
		rast = raster(as.matrix(c(min(decimalDates),max(decimalDates)))); par(lwd=0.2, col="gray30")
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.500,0.512,0.40,0.750),
			 legend.args=list(text="", cex=0.8, line=0.3, col="gray30"), horizontal=F,
		     axis.args=list(cex.axis=0.8, lwd=0, lwd.tick=0.2, tck=-0.6, col.axis="gray30", line=0, mgp=c(0,0.5,0)))
	}
cols = colourScale[colIndices]; metadata[,4] = decimalDates; metadata[,5] = cols
if (writingFiles) write.table(metadata, "LSVD_all_alignment4.csv", col.names=F, row.names=F, quote=F, sep=";")
metadata = metadata[,1:2]
for (i in 1:dim(metadata)[1])
	{
		metadata[i,1] = substr(metadata[i,1],1,6)
	}
if (writingFiles) write.table(metadata, "LSVD_Network_pops.csv", col.names=F, row.names=F, quote=F, sep=";")
populations = cbind(different_countries, different_colours)
for (i in 1:dim(populations)[1])
	{
		populations[i,1] = substr(populations[i,1],1,2)
	}
if (writingFiles) write.table(populations, "LSVD_Network_cols.csv", col.names=F, row.names=F, quote=F, sep=";")

# 2. Testing for a signal of recombination within the global alignment

	# --> the phi-test did find statistically significant evidence for recombination (p < 0.001)

# 3. Generating an haplotype network for the global data set (Network)

source("networkGraph_SPADS.r")
network_out = "LSVD_all_alignment4.out"
populations = "LSVD_Network_pops.csv"
colour_code = "LSVD_Network_cols.csv"
scaleFactor1 = 4; scaleFactor2 = 1
networkGraph_SPADS(network_out, populations, colour_code, scaleFactor1, scaleFactor2)

# 4. Inferring a maximum likelihood phylogeny for the global data set

system("IQTREE_1.6.12_on_MacOS/iqtree -s LSVD_all_alignment1.fas -m MFP -mem 10Go -mset GTR -nt 4 -b 100")

tre = readAnnotatedNexus("LSVD_all_alignment1.tre"); tre$tip.label = paste0("  ",gsub("'","",tre$tip.label))
tab = read.csv("LSVD_all_alignment1.csv", head=T, sep=";"); tab[which(tab[,"country"]=="China "),"country"] = "China"
for (i in 1:length(tre$tip.label))
	{
		country = tab[which(grepl(gsub("  ","",tre$tip.label[i]),tab[,"seqID"])),"country"]
		tre$tip.label[i] = paste0(tre$tip.label[i]," (",country,")")
	}
dev.new(width=9, height=6); par(mar=c(0.7,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o"); indices = c()
plot(tre, show.tip.label=T, show.node.label=F, edge.width=0.75, cex=0.5, col="gray30", edge.color="gray30")
tre = readAnnotatedNexus("LSVD_all_alignment1.tre"); tre$tip.label = paste0("  ",gsub("'","",tre$tip.label))
for (i in 1:dim(tre$edge)[1])
	{
		if (!tre$edge[i,2]%in%tre$edge[,1])
			{
				index = which(grepl(gsub("  ","",tre$tip.label[tre$edge[i,2]]),tab[,"seqID"]))
				colour = different_colours[which(different_countries==tab[index,"country"])]
				nodelabels(node=tre$edge[i,2], pch=16, cex=0.70, col=colour)
				nodelabels(node=tre$edge[i,2], pch=1, cex=0.70, col="gray30", lwd=0.2)
				# tiplabels(tab[index,"seqID"], tip=tre$edge[i,2], col="gray30", cex=0.50, frame="none", adj=-0.2)
			}	else	{
				if (tre$annotations[[i]]$label >= 70)
					{
						nodelabels(node=tre$edge[i,2], pch=16, cex=0.50, col="gray30")
					}
			}
	}
add.scale.bar(x=0.0, y=-0.1, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)

# 5. Generating an overall sampling map for all the samples

tab = read.csv("LSVD_all_alignment1.csv", head=T, sep=";"); tab[which(tab[,"country"]=="China "),"country"] = "China"
unique_coordinates = unique(tab[,"latitude_longitude"]); e_Palearctic = extent(-25, 180, -40, 75)
countries = crop(gBuffer(shapefile("World_countries_shapefile/World_countries_shapefile.shp") , byid=T, width=0), e_Palearctic)
borders = crop(shapefile("Only_international_borders/Only_international_borders.shp"), e_Palearctic)
coasts = crop(shapefile("Only_coast_lines_borders/Only_coast_lines_borders.shp"), e_Palearctic)
sampling = as.data.frame(matrix(nrow=length(unique_coordinates), ncol=5))
colnames(sampling) = c("longitude","latitude","country","colour","samples")
for (i in 1:dim(sampling)[1])
	{
		sampling[i,1] = as.numeric(unlist(strsplit(unique_coordinates[i],", "))[2])
		sampling[i,2] = as.numeric(unlist(strsplit(unique_coordinates[i],", "))[1])
		indices = which(tab[,"latitude_longitude"]==unique_coordinates[i])
		sampling[i,3] = tab[indices[1],"country"]; sampling[i,5] = length(indices)
		index = which(different_countries==tab[indices[1],"country"])
		sampling[i,4] = different_colours[index]
	}
sampling = sampling[which(!is.na(sampling[,1])),]

pdf("LSDV_sampling_map.pdf", width=8, height=5)
par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
plot(countries, col="gray90", border=NA, ann=F, axes=F)
plot(borders, col="white", lwd=0.3, add=T)
plot(coasts, col="gray70", lwd=0.5, add=T)
for (i in 1:dim(sampling)[1])
	{
		cex = (((sqrt(sampling[i,5]/pi))*2)/((sqrt(1/pi))*2))*0.5
		points(sampling[i,1:2], col=sampling[i,4], cex=cex, pch=16)
		points(sampling[i,1:2], cex=cex, pch=1, lwd=0.5, col="gray30")
	}
rect(-25, -40, 180, 75, lwd=0.2, border="gray30")
dev.off()

# 6. Inferring the recombination breakpoints with GARD (www.datamonkey.org/GARD)

	# - parameters: site-to-site rate variation = beta-gamma; rate classes = 6
	# - results: http://www.datamonkey.org/gard/6424343907f53a24047dbbbf
	# --> NRRs: 1-189, 190-1489, 1490-1660, 1661-1827, 1828-1955, 1956-2088

# 7. Preparing the input files for the analyses to perform in SPADS

selected_samples = gsub("'","",read.nexus("LSVD_all_1_wild_type.tre")$tip.label)
R4_clade_samples = gsub("'","",read.nexus("LSVD_all_1_R4_clade.tre")$tip.label)
african_clade = c("OP688128.1","ON400507.1","OP297402.1","OK422492.1","OP688129.1","KX683219.1","AF325528.1","MN072619.1")
txt = scan(paste0("LSVD_all_alignment2.phy"), what="", sep="\n", quiet=T)
txt = gsub("6-16-Gr","6_16-Gr",txt); IDs = c(); indices1 = c(); indices2 = c()
tab = read.csv("LSVD_all_alignment1.csv", head=T, sep=";"); indices3 = c()
for (i in 2:length(txt)) IDs = c(IDs, unlist(strsplit(txt[i]," "))[1])
for (i in 1:length(IDs))
	{
		ID1 = gsub("6_16","6-16",unlist(strsplit(IDs[i],"-"))[1])
		ID2 = unlist(strsplit(IDs[i],"-"))[length(unlist(strsplit(IDs[i],"-")))]
		index = which((grepl(ID1,tab[,"seqID"]))&(grepl(ID2,tab[,"country"])))
		if (length(index) != 1) cat("Problem to find sample ",IDs[i]," in metadata\n",sep="")
		ID3 = unlist(strsplit(tab[index,"seqID"]," "))[1]
		if ((ID3%in%selected_samples)&(!ID3%in%african_clade)) indices1 = c(indices1, i+1)
		if (ID3%in%selected_samples) indices2 = c(indices2, i+1)
		if (ID3%in%R4_clade_samples) indices3 = c(indices3, i+1)
	}
txt1 = txt[c(1,indices1)]; txt1[1] = paste(length(txt1)-1,unlist(strsplit(txt1[1]," "))[2],"1")
txt2 = txt[c(1,indices2)]; txt2[1] = paste(length(txt2)-1,unlist(strsplit(txt2[1]," "))[2],"1")
txt3 = txt[c(1,indices3)]; txt3[1] = paste(length(txt3)-1,unlist(strsplit(txt3[1]," "))[2],"1")
if (writingFiles)
	{
		write(txt1, "LSVD_all_SPADS_analyses/Analysis_1/LSVD_all_alignment2_mod.phy")
		write(txt2, "LSVD_all_SPADS_analyses/Analysis_2/LSVD_all_alignment2_mod.phy")
		write(txt3, "LSVD_all_SPADS_analyses/Analysis_3/LSVD_all_alignment2_mod.phy")
	}
for (h in 1:3)
	{
		txt = scan(paste0("LSVD_all_SPADS_analyses/Analysis_",h,"/LSVD_all_alignment2_mod.phy"), what="", sep="\n", quiet=T)
		txt = gsub("6-16-Gr","6_16-Gr",txt); IDs = c(); pops = c()
		tab = read.csv("LSVD_all_alignment1.csv", head=T, sep=";")
		for (i in 2:length(txt)) IDs = c(IDs, unlist(strsplit(txt[i]," "))[1])
		coordinates = matrix(nrow=length(IDs), ncol=4)
		for (i in 1:length(IDs))
			{
				ID1 = gsub("6_16","6-16",unlist(strsplit(IDs[i],"-"))[1])
				ID2 = unlist(strsplit(IDs[i],"-"))[length(unlist(strsplit(IDs[i],"-")))]
				index = which((grepl(ID1,tab[,"seqID"]))&(grepl(ID2,tab[,"country"])))
				if (length(index) != 1) cat("Problem to find sample ",IDs[i]," in metadata\n",sep="")
				longitude = as.numeric(unlist(strsplit(tab[index,"latitude_longitude"],", "))[2])
				latitude = as.numeric(unlist(strsplit(tab[index,"latitude_longitude"],", "))[1])
				coordinates[i,] = cbind(tab[index,"latitude_longitude"], longitude, latitude, tab[index,"country"])
			}
		unique_coordinates = unique(coordinates[,1]); old_IDs = IDs
		for (i in 1:length(IDs))
			{
				ID1 = unlist(strsplit(IDs[i],"-"))[1]
				ID2 = which(unique_coordinates==coordinates[i,1])
				if (nchar(ID2) == 1) ID2 = paste0("0",ID2)
				IDs[i] = paste0(ID1,"-",ID2)
			}
		for (i in 1:length(IDs)) pops = c(pops, unlist(strsplit(IDs[i],"-"))[2])
		pops = unique(pops); populations = matrix(nrow=length(pops), ncol=3)
		groups = matrix(nrow=length(pops), ncol=2)
		for (i in 1:length(pops))
			{
				populations[i,1] = pops[i]
				populations[i,2:3] = coordinates[as.numeric(pops[i]),2:3]
				groups[i,] = cbind(pops[i], pops[i])
			}
		for (i in 1:length(old_IDs))
			{
				txt = gsub(old_IDs[i], IDs[i], txt)
			}
		indices = which(is.na(populations[,2]))
		if (length(indices) > 0)
			{
				populations_to_discard = populations[indices,1]
				populations = populations[-indices,]
				indices = which(groups[,1]==populations_to_discard)
				groups = groups[-indices,]; indices = c()
				for (i in 1:length(populations_to_discard))
					{
						indices = c(indices, which(grepl(paste0("-",populations_to_discard[i]," "),txt)))
					}
				txt = txt[-indices]
			}
		if (writingFiles)
			{
				write(txt, paste0("LSVD_all_SPADS_analyses/Analysis_",h,"/LSVD_all_alignment2_locus1.phy"))
				write.table(populations, paste0("LSVD_all_SPADS_analyses/Analysis_",h,"/LSVD_all_alignment2_populations.txt"), sep=" ", row.names=F, col.names=F, quote=F)
				write.table(groups, paste0("LSVD_all_SPADS_analyses/Analysis_",h,"/LSVD_all_alignment2_groups.txt"), sep=" ", row.names=F, col.names=F, quote=F)		
			}
	}
	
# 8. Testing for a signal of recombination within the two wild type clades

	# --> the phi-test did not find statistically significant evidence for recombination (p = 0.644)

# 9. Investigating the isolation-by-distance pattern based on the IID2 metric

	# 9.1. When considering all the two wild type clades

genDis = as.matrix(read.table("LSVD_all_SPADS_analyses/Analysis_1/LSVD_all_alignment2_GDisPAL_input_distances_matrix_IID2.txt"))
coordinates = read.table("LSVD_all_SPADS_analyses/Analysis_1/LSVD_all_alignment2_GDisPAL_input_coordinates.txt")[,2:1]
geoDis = matrix(0, nrow=dim(coordinates)[1], ncol=dim(coordinates)[1])
for (i in 2:dim(geoDis)[1])
	{
		x1 = cbind(coordinates[i,1],coordinates[i,2])
		for (j in 1:(i-1))
			{
				x2 = cbind(coordinates[j,1],coordinates[j,2])
				d = rdist.earth(x1, x2, miles=F, R=NULL)
				geoDis[i,j] = log(d+1); geoDis[j,i] = log(d+1)
			}
	}
ape::mantel.test(geoDis, genDis, nperm=999) # Z-stat = 0.211, p-value = 0.001
vegan::mantel(geoDis, genDis, method="spearman", permutations=999) # rS = 0.648, p-value = 0.001
geoDis = geoDis[lower.tri(geoDis)]; genDis = genDis[lower.tri(genDis)]
lr = lm(genDis ~ geoDis); summary(lr) # R2 = 0.37 (p-value << 0.001)
shapiro.test(residuals(lr)) # plot(lr, which=2) # p-value << 0.001

	# 9.2. When considering only the main wild type clade

genDis = as.matrix(read.table("LSVD_all_SPADS_analyses/Analysis_2/LSVD_all_alignment2_GDisPAL_input_distances_matrix_IID2.txt"))
coordinates = read.table("LSVD_all_SPADS_analyses/Analysis_2/LSVD_all_alignment2_GDisPAL_input_coordinates.txt")[,2:1]
geoDis = matrix(0, nrow=dim(coordinates)[1], ncol=dim(coordinates)[1])
for (i in 2:dim(geoDis)[1])
	{
		x1 = cbind(coordinates[i,1],coordinates[i,2])
		for (j in 1:(i-1))
			{
				x2 = cbind(coordinates[j,1],coordinates[j,2])
				d = rdist.earth(x1, x2, miles=F, R=NULL)
				geoDis[i,j] = log(d+1); geoDis[j,i] = log(d+1)
			}
	}
ape::mantel.test(geoDis, genDis, nperm=999) # Z-stat = 1.157, p-value = 0.001
vegan::mantel(geoDis, genDis, method="spearman", permutations=999) # rS = 0.762, p-value = 0.001
geoDis = geoDis[lower.tri(geoDis)]; genDis = genDis[lower.tri(genDis)]
lr = lm(genDis ~ geoDis); summary(lr) # R2 = 0.39 (p-value << 0.001)
shapiro.test(residuals(lr)) # plot(lr, which=2) # p-value << 0.001

	# 9.3. When excluding African, Kazachstan and Russian samples

toKeep = which((coordinates[,1]<45)&(coordinates[,2]>25)); geoDis_sub = geoDis[toKeep,toKeep]; genDis_sub = genDis[toKeep,toKeep]
ape::mantel.test(geoDis_sub, genDis_sub, nperm=9999) # Z-stat = 0.088, p-value = 0.001
vegan::mantel(geoDis_sub, genDis_sub, method="spearman", permutations=9999) # rS = 0.515, p-value = 0.002
geoDis_sub = geoDis_sub[lower.tri(geoDis_sub)]; genDis_sub = genDis_sub[lower.tri(genDis_sub)]
lr = lm(genDis_sub ~ geoDis_sub); summary(lr) # R2 = 0.42 (p-value << 0.001)
shapiro.test(residuals(lr)) # plot(lr, which=2) # p-value << 0.001

	# 9.4. When only considering the European sequences (minus Turkey)

toKeep = which((coordinates[,1]<30)&(coordinates[,2]>30)); geoDis_sub = geoDis[toKeep,toKeep]; genDis_sub = genDis[toKeep,toKeep]
ape::mantel.test(geoDis_sub, genDis_sub, nperm=9999) # Z-stat = 0.019, p-value = 0.266
vegan::mantel(geoDis_sub, genDis_sub, method="spearman", permutations=9999) # rS = 0.100, p-value = 0.239
geoDis_sub = geoDis_sub[lower.tri(geoDis_sub)]; genDis_sub = genDis_sub[lower.tri(genDis_sub)]
lr = lm(genDis_sub ~ geoDis_sub); summary(lr) # R2 = 0.01 (p-value = 0.164)
shapiro.test(residuals(lr)) # plot(lr, which=2) # p-value < 0.001

	# 9.5. When considering the "R4" clade (see the Network/tree figure)

genDis = as.matrix(read.table("LSVD_all_SPADS_analyses/Analysis_3/LSVD_all_alignment2_GDisPAL_input_distances_matrix_IID2.txt"))
coordinates = read.table("LSVD_all_SPADS_analyses/Analysis_3/LSVD_all_alignment2_GDisPAL_input_coordinates.txt")[,2:1]
geoDis = matrix(0, nrow=dim(coordinates)[1], ncol=dim(coordinates)[1])
for (i in 2:dim(geoDis)[1])
	{
		x1 = cbind(coordinates[i,1],coordinates[i,2])
		for (j in 1:(i-1))
			{
				x2 = cbind(coordinates[j,1],coordinates[j,2])
				d = rdist.earth(x1, x2, miles=F, R=NULL)
				geoDis[i,j] = log(d+1); geoDis[j,i] = log(d+1)
			}
	}
ape::mantel.test(geoDis, genDis, nperm=999) # Z-stat = 0.094, p-value = 0.721
vegan::mantel(geoDis, genDis, method="spearman", permutations=999) # rS = 0.089, p-value = 0.270
geoDis = geoDis[lower.tri(geoDis)]; genDis = genDis[lower.tri(genDis)]
lr = lm(genDis ~ geoDis); summary(lr) # R2 = -0.008 (p-value = 0.826)
shapiro.test(residuals(lr)) # plot(lr, which=2) # p-value << 0.001

if (savingPlots)
	{
		pdf("Isolation-by-distance.pdf", width=3.5, height=3.7) # dev.new(width=3.5, height=3.7)
		par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(3.3,3.5,2,2), lwd=0.2, col="gray30")
		plot(geoDis, genDis, col=NA, pch=16, cex=0.8, axes=F, ann=F, frame=T)
		abline(lr, lwd=1, col="gray50")
		for (i in 1:length(geoDis))
			{
				points(geoDis[i], genDis[i], col="gray90", pch=16, cex=0.8)
				points(geoDis[i], genDis[i], col="gray30", pch=1, cex=0.8, lwd=0.2)
			}
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.0, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.30,0), lwd=0.0, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30")
		title(ylab="genetic distance (IID2)", cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
		title(xlab="geographic distance (km, log-transformed)", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		dev.off()
	}

# 10. Investigating the temporal signal and preparing the BEAST input files

tre = read.nexus("LSVD_all_1_wild_type.tre")
tre$tip.label = gsub("'","",tre$tip.label)
tab = read.table("LSVD_all_alignment1.csv", head=T, sep=";")
mat = matrix(nrow=length(tre$tip.label), ncol=5)
colnames(mat) = c("trait","date","latitude","longitude","network_ID")
mat[,"trait"] = tre$tip.label
for (i in 1:dim(mat)[1])
	{
		index = which(grepl(tre$tip.label[i],tab[,"seqID"]))
		txt = gsub("May","05",gsub("\\/","-",tab[index,"sampling_date"]))
		coordinates = tab[index,"latitude_longitude"]
		if (length(unlist(strsplit(txt,"-"))) == 1)
			{
				date = txt
			}
		if (length(unlist(strsplit(txt,"-"))) == 2)
			{
				year = unlist(strsplit(txt,"-"))[2]
				month = unlist(strsplit(txt,"-"))[1]
				date = paste(year,month,sep="-")
			}
		if (length(unlist(strsplit(txt,"-"))) == 3)
			{
				year = unlist(strsplit(txt,"-"))[3]
				month = unlist(strsplit(txt,"-"))[2]
				day = unlist(strsplit(txt,"-"))[1]
				date = paste(year,month,day,sep="-")
			}		
		mat[i,"date"] = date
		mat[i,"latitude"] = unlist(strsplit(coordinates,", "))[1]
		mat[i,"longitude"] = unlist(strsplit(coordinates,", "))[2]
		country = substr(tab[index,"country"],1,2)
		ID = unlist(strsplit(tab[index,"seqID"],"\\."))[1]
		ID = substr(ID,nchar(ID)-3,nchar(ID))
		mat[i,"network_ID"] = paste0(ID,"-",country)
	}
write.table(mat, "LSVD_all_1_wild_type.txt", row.names=F, quote=F, sep="\t")

mat = mat[which((!is.na(mat[,"latitude"]))&(!is.na(mat[,"longitude"]))),]
txt = scan(paste0("LSVD_all_alignment2.phy"), what="", sep="\n", quiet=T)
txt = txt[2:length(txt)]; buffer = c()
for (i in 1:length(txt))
	{
		network_ID = unlist(strsplit(txt[i],"   "))[1]
		sequence = unlist(strsplit(txt[i],"   "))[2]
		index = which(mat[,"network_ID"]==network_ID)
		if (length(index) == 1)
			{
				fasta_ID = mat[index,"trait"]
				buffer = c(buffer,paste0(">",fasta_ID),sequence)
			}
	}
write.table(mat, "LSVD_all_1_WTs_BEAST/LSVD_all_1_wild_type.txt", row.names=F, quote=F, sep="\t")
write(buffer, "LSVD_all_1_WTs_BEAST/LSVD_all_1_wild_type.fas")

files = list.files("LSVD_all_1_WTs_BEAST/LSVD_all_1_WTs_polygons/Longitude_latitude/")
for (i in 1:length(files))
	{
		txt = scan(paste0("LSVD_all_1_WTs_BEAST/LSVD_all_1_WTs_polygons/Longitude_latitude/",files[i]), what="", sep="\n", quiet=T); buffer = txt
		for (j in 1:length(txt))
			{
				if (grepl("\t\t\t",txt[j]))
					{
						line = unlist(strsplit(gsub("\t\t\t","",txt[j]),","))
						buffer[j] = paste0("\t\t\t",line[2],",",line[1],",",line[3])
					}
			}
		write(buffer, paste0("LSVD_all_1_WTs_BEAST/LSVD_all_1_WTs_polygons/",files[i]))
	}

files = list.files("LSVD_all_1_WTs_BEAST/LSVD_all_1_WTs_polygons/"); samples = gsub(".kml","",files)
for (i in 1:dim(mat)[1])
	{
		if (mat[i,"trait"]%in%samples)
			{
				index = which(samples==mat[i,"trait"])
				txt = scan(paste0("LSVD_all_1_WTs_BEAST/LSVD_all_1_WTs_polygons/",files[index]), what="", sep="\n", quiet=T)
				txt = gsub("\t\t\t","",txt[which(grepl("\t\t\t",txt))])
				tmp = matrix(nrow=length(txt), ncol=2)
				for (j in 1:length(txt))
					{
						tmp[j,1] = as.numeric(unlist(strsplit(txt[j],","))[2]); tmp[j,2] = as.numeric(unlist(strsplit(txt[j],","))[1])
					}
				x = as.numeric(mat[i,"longitude"]); y = as.numeric(mat[i,"latitude"])
				if (point.in.polygon(x,y,tmp[,1],tmp[,2]) != 1) print(c(i,mat[i,"trait"]))
			}
	}

directory = "LSVD_all_1_WTs_polygons"
files = list.files("LSVD_all_1_WTs_BEAST/LSVD_all_1_WTs_polygons/")
files = files[which(grepl(".kml",files))]; samples = gsub(".kml","",files)
xml = scan("LSVD_all_1_WTs_BEAST/LSVD_all_1_wild_type.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
xml = gsub("LSVD_all_1_wild_type","LSVD_all_1_WTs_pols",xml)
sink(file=paste0("LSVD_all_1_WTs_BEAST/LSVD_all_1_WTs_pols.xml"))
for (i in 1:length(xml))
	{
		cat(xml[i]); cat("\n")
		if (xml[i] == "\t<!-- INSERT leafTraitParameter -->")
			{
				cat("\n")
				for (j in 1:length(samples))
					{
						cat(paste("\t<leafTraitParameter id=\"",samples[j],".trait\" taxon=\"",samples[j],"\">",sep="")); cat("\n")
						cat(paste("\t\t<treeModel idref=\"treeModel\"/>",sep="")); cat("\n")
						cat(paste("\t\t<parameter idref=\"leaf.location\"/>",sep="")); cat("\n")
						cat(paste("\t</leafTraitParameter>",sep="")); cat("\n")
					}
				cat("\n")
				for (j in 1:length(samples))
					{
						cat(paste("\t<flatGeoSpatialPrior id=\"",samples[j],"_polygon\" taxon=\"",samples[j],"\" kmlFileName=\"",directory,"/",samples[j],".kml\" inside=\"true\" union=\"true\" cache=\"true\">",sep="")); cat("\n")
						cat(paste("\t\t<data>",sep="")); cat("\n")
						cat(paste("\t\t\t<parameter idref=\"",samples[j],".trait\"/>",sep="")); cat("\n")
						cat(paste("\t\t</data>",sep="")); cat("\n")
						cat(paste("\t</flatGeoSpatialPrior>",sep="")); cat("\n")
					}
				cat("\n")		
			}
		if (xml[i]=="\t\t<!-- INSERT uniformGeoSpatialOperator -->")
			{
				cat("\n")
				for (j in 1:length(samples))
					{
						cat(paste("\t\t<uniformGeoSpatialOperator weight=\"0.01\">",sep="")); cat("\n")
						cat(paste("\t\t\t<parameter idref=\"",samples[j],".trait\"/>",sep="")); cat("\n")
						cat(paste("\t\t\t<flatGeoSpatialPrior idref=\"",samples[j],"_polygon\"/>",sep="")); cat("\n")
						cat(paste("\t\t</uniformGeoSpatialOperator>",sep="")); cat("\n")
					}
				cat("\n")
			}
		if (xml[i]=="\t\t\t\t<!-- INSERT geoDistributionCollection -->")
			{
				cat("\n")
				cat("\t\t\t\t<geoDistributionCollection id=\"allGeoDistributions\">"); cat("\n")
				for (j in 1:length(samples))
					{
						cat(paste("\t\t\t\t<flatGeoSpatialPrior idref=\"",samples[j],"_polygon\"/>",sep="")); cat("\n")
					}
				cat("\t\t\t\t</geoDistributionCollection>"); cat("\n")
				cat("\n")
			}					
	}
sink(NULL)

	# Notes:
		# - to avoid having to use a small jitter, KY829023 was "commented" in the XML files because coming from the exact same municipality as GRC478
		# - for the same reason, the "starting" coordinates assigned to the sequences from Bangladesh and Kenya in the XML files have been manually re-sampled
		
# 11. Extracting the spatio-temporal information embedded in posterior trees

analyses = c("LSVD_all_1_WTs_pols")
nberOfTreesToSample = 1000; burnIn = 1001
for (i in 1:length(analyses))
	{
		log = scan(paste0("LSVD_all_1_WTs_BEAST/",analyses[i],".log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		index1 = 6+burnIn; index2 = length(log); interval = round((index2-index1)/nberOfTreesToSample)
		indices = seq(index2-((nberOfTreesToSample-1)*interval),index2,interval)
		write(log[c(5,indices)], paste0(analyses[i],".log"))
		trees = scan(paste0("LSVD_all_1_WTs_BEAST/",analyses[i],".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		index1 = which(trees=="\t\t;")[length(which(trees=="\t\t;"))]
		index2 = index1 + burnIn + 1
		indices3 = which(grepl("tree STATE",trees)); index3 = indices3[length(indices3)]
		interval = floor((index3-(index1+burnIn))/nberOfTreesToSample)
		indices = seq(index3-((nberOfTreesToSample-1)*interval),index3,interval)
		selected_trees = c(trees[c(1:index1,indices)],"End;")
		write(selected_trees, paste0(analyses[i],".trees"))
	}
		# To do: getting and annotating the MCC tree with TreeAnnotator, using the 1,000 selected trees as an input
		# Tracer: ucld.mean = 5.91E-5, 95%HPD = [5.41E-7, 3.12E-4]; root age = 1888.1, 95% HPD = [1752.9-1957.4]

source("Tree_data_extraction1.r") # for the MCC tree
source("Tree_data_extraction2.r") # for the posterior trees
mostRecentSamplingDates = rep(NA, length(analyses))
mostRecentSamplingDatum = 2021.5 # ?? (unprecise...)
for (i in 1:length(analyses))
	{
		mcc_tre = readAnnotatedNexus(paste0(analyses[i],".tree"))
		mostRecentSamplingDates[i] = mostRecentSamplingDatum
		mcc_tab = Tree_data_extraction1(mcc_tre, mostRecentSamplingDatum)
		write.csv(mcc_tab, paste0(analyses[i],".csv"), row.names=F, quote=F)
	}
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0(analyses[i],"_ext"); mostRecentSamplingDatum = mostRecentSamplingDates[i]
		allTrees = readAnnotatedNexus(paste0(analyses[i],".trees"))
		for (j in 1:length(allTrees))
			{
				tab = Tree_data_extraction2(allTrees[[j]], mostRecentSamplingDatum)
				write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
			}
	}

# 12. Mapping the dispersal history of LSVD lineages of the wild type clades

analysis = c("LSVD_all_1_WTs_pols"); nberOfExtractionFiles = 1000
mostRecentSamplingDatum = 2021.5 # ?? (unprecise...)
localTreesDirectory = paste0(analysis,"_ext")
rootHeights = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
	{
		csv = read.csv(paste(localTreesDirectory,"/TreeExtractions_",i,".csv",sep=""), header=T)
		rootHeights[i] = min(csv[,"startYear"])
	}
mcc = read.csv(paste0(analysis,".csv"), head=T)
mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
minYear = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]; maxYear = mostRecentSamplingDatum
colour_scale = met.brewer(name="Hiroshige", n=111, type="continuous")[1:101]

tree = readAnnotatedNexus(paste0(analysis,".tree")); tree$tip.label = gsub("'","",tree$tip.label)
rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingDatum-rootHeight
minYear = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]; maxYear = mostRecentSamplingDatum

pdf(paste0(analysis,"_NEW1.pdf"), width=3.5, height=4.05); # dev.new(width=3.5, height=4.05)
par(mar=c(1,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o", col="gray30"); plottingRootNode = TRUE
plot(tree, show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, 
	 x.lim=c(minYear-(maxYear-max(nodeHeights(tree))), max(nodeHeights(tree))), col="gray30", edge.color="gray30")
tree_obj = get("last_plot.phylo", envir=.PlotPhyloEnv); plottingRootBar = FALSE
for (j in 1:dim(tree$edge)[1])
	{
		endYear = root_time+nodeHeights(tree)[j,2]
		endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
		endYear_colour = colour_scale[endYear_index]
		if ((tree$edge[j,2]%in%tree$edge[,1])&(length(tree$annotations[[j]]$`height_95%_HPD`) > 1))
			{
				x1 = (mostRecentSamplingDatum-tree$annotations[[j]]$`height_95%_HPD`[[2]])-root_time
				x2 = (mostRecentSamplingDatum-tree$annotations[[j]]$`height_95%_HPD`[[1]])-root_time
				lines(x=c(x1,x2), y=rep(tree_obj$yy[tree$edge[j,2]],2), lwd=3.5, lend=0, col=paste0(endYear_colour,"40"))
			}
		if ((plottingRootBar == FALSE)&&(!tree$edge[j,1]%in%tree$edge[,2]))
			{
				endYear = root_time+nodeHeights(tree)[j,1]
				endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				endYear_colour = colour_scale[endYear_index]
				x1 = (mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]])-root_time
				x2 = (mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[1]])-root_time
				lines(x=c(x1,x2), y=rep(tree_obj$yy[tree$edge[j,1]],2), lwd=3.5, lend=0, col=paste0(endYear_colour,"40"))
				plottingRootBar = TRUE
			}				
	}
for (j in 1:dim(tree$edge)[1])
	{
		endYear = root_time+nodeHeights(tree)[j,2]
		endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
		endYear_colour = colour_scale[endYear_index]
		if ((tree$edge[j,2]%in%tree$edge[,1])&&(tree$annotations[[j]]$posterior >= 0.95))
			{
				nodelabels(node=tree$edge[j,2], pch=16, cex=0.50, col=endYear_colour)
				nodelabels(node=tree$edge[j,2], pch=1, cex=0.50, col="gray30", lwd=0.2)
			}
		if (!tree$edge[j,2]%in%tree$edge[,1])
			{
				nodelabels(node=tree$edge[j,2], pch=16, cex=0.50, col=endYear_colour)
				nodelabels(node=tree$edge[j,2], pch=1, cex=0.50, col="gray30", lwd=0.2)
			}
		if ((plottingRootNode == TRUE)&&(!tree$edge[j,1]%in%tree$edge[,2]))
			{
				endYear = root_time+nodeHeights(tree)[j,1]
				endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				endYear_colour = colour_scale[endYear_index]; plottingRootNode = TRUE		
				nodelabels(node=tree$edge[j,1], pch=16, cex=0.50, col=endYear_colour)
				nodelabels(node=tree$edge[j,1], pch=1, cex=0.50, col="gray30", lwd=0.2)
			}
	}
selectedDates = seq(1750,2000,50); selectedLabels = selectedDates
selectedDates = c(minYear, selectedDates, maxYear); selectedLabels = c("", selectedLabels, "")
axis(lwd=0.3, at=selectedDates-root_time, labels=selectedLabels, cex.axis=0.60, mgp=c(0,-0.15,-0.3), lwd.tick=0.3, 
	 col.lab="gray30", col="gray30", tck=-0.013, side=1)
dev.off()

e_Palearctic = extent(-25, 125, -40, 75)
countries = crop(gBuffer(shapefile("World_countries_shapefile/World_countries_shapefile.shp") , byid=T, width=0), e_Palearctic)
borders = crop(shapefile("Only_international_borders/Only_international_borders.shp"), e_Palearctic)
coasts = crop(shapefile("Only_coast_lines_borders/Only_coast_lines_borders.shp"), e_Palearctic)

prob = 0.80; precision = 25; startDatum=minYear
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))

endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
		polygons_colours[i] = paste0(colour_scale[polygon_index],"20")
	}

cutOffs = c(maxYear); croppingPolygons = FALSE
for (h in 1:length(cutOffs))
	{
		pdf(paste0(analysis,"_NEW2.pdf"), width=5, height=4.05)
		par(mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o")
		plot(countries, col="gray90", border=NA, ann=F, axes=F)
		plot(borders, col="white", lwd=0.3, add=T)
		plot(coasts, col="gray70", lwd=0.5, add=T)
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
		if (h == length(cutOffs))
			{
				plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.50,0.85,0.11,0.12),
	 				 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
					 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-1.0, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.15,0)))
			 }
		for (i in 1:length(polygons))
			{
				if (as.numeric(names(polygons[[i]])) <= cutOffs[h])
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(countries)
						if (croppingPolygons == TRUE) pol = crop(pol, countries)
						plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
					}
			}
		for (i in 1:dim(mcc)[1])
			{
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				  		  		    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
				  	}
			}
		for (i in dim(mcc)[1]:1)
			{
				if ((mcc[i,"startYear"] <= cutOffs[h])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
					{
						startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.4)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
					}
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.4)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
					}
			}
		rect(-25, -40, 125, 75, lwd=0.2, border="gray30")
		dev.off()
	}

