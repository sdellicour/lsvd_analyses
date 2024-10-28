networkGraph_SPADS = function(network_out, populations, colour_code, scaleFactor1, scaleFactor2) {

	network_txt = scan(network_out, what="", sep="\n", quiet=T)
	partition = read.csv(populations, sep=";", header=F, colClasses="character", stringsAsFactors=F)
	colours = read.csv(colour_code, sep=";", header=F, colClasses="character", stringsAsFactors=F)
	for (i in 1:dim(partition)[1])
		{
			partition[i,2] = colours[which(colours[,1]==partition[i,2]),2]
		}
	indices = which(grepl("Link ",network_txt))
	links = matrix(nrow=length(indices), ncol=4)
	colnames(links) = c("source","target","mutations","label")
	nodeIDs = c()
	for (i in 1:length(indices))
		{
			line = network_txt[indices[i]]
			line = gsub("mv +","mv",line)
			line = gsub(",","",unlist(strsplit(line," +")))
			index1 = which(line=="Seq.")[1]+1
			index2 = which(line=="Seq.")[2]+1
			index3 = which(line=="at")
			node1 = line[index1]
			node2 = line[index2]
			nodeIDs = c(nodeIDs, node1, node2)
			mutations = (length(line)-index3)/2
			label = mutations
			if (label == 1) label = NA
			links[i,] = cbind(node1,node2,mutations,label)
		}
	nodeIDs = unique(nodeIDs); nberOfNodes = length(nodeIDs)
	nodes = matrix(nrow=nberOfNodes, ncol=2)
	indices = which(grepl("Tax.",network_txt))
	grouping = matrix(nrow=length(indices), ncol=2)
	colnames(nodes) = c("sequences","diameter")
	row.names(nodes) = nodeIDs
	nodes[,"sequences"] = rep(1, dim(nodes)[1])
	for (i in 1:length(indices))
		{
			line = network_txt[indices[i]]
			line = gsub(",","",unlist(strsplit(line," +")))
			index1 = which(grepl("ax.",line))[1]+1
			index2 = which(grepl("ax.",line))[2]+1
			nodes[line[index1],"sequences"] = nodes[line[index1],"sequences"]+1
			grouping[i,] = cbind(line[index1],line[index2])
		}
	for (i in 1:dim(nodes)[1])
		{
			area = nodes[i,"sequences"]
			nodes[i,"diameter"] = sqrt(area/pi)*2*scaleFactor1
		}
	nodes[which(grepl("mv",row.names(nodes))),"diameter"] = 0
	big_nodes = unique(grouping[,1])
	nodes_list1 = list()
	for (i in 1:length(big_nodes))
		{
			nodes_list1[[i]] = c(big_nodes[i], grouping[which(grouping[,1]==big_nodes[i]),2])
		}
	small_nodes = nodeIDs[!nodeIDs%in%big_nodes]
	for (i in 1:length(small_nodes))
		{
			nodes_list1[[length(big_nodes)+i]] = small_nodes[i]
		}
	indices = which(grepl("Mapping of taxon; ",network_txt))
	nodeIDs_correspondence = matrix(nrow=length(indices), ncol=2)
	for (i in 1:length(indices))
		{
			line = gsub("mv +","mv",network_txt[indices[i]])
			line = unlist(strsplit(line," +"))
			line = unlist(strsplit(line[length(line)],";"))
			nodeIDs_correspondence[i,] = cbind(line[1],line[2])
		}
	nodes_list2 = list()
	for (i in 1:dim(nodes)[1])
		{
			index = c()
			for (j in 1:length(nodes_list1))
				{
					if (nodes_list1[[j]][1] == row.names(nodes)[i]) index = c(index, j)
				}
			if (length(index) != 1)
				{
					print(i)
				}	else	{
					nodes_list2[[i]] = nodes_list1[[index]]
				}
		}
	nodes_values = list(); nodes_cols = list()
	for (i in 1:length(nodes_list2))
		{
			seqs1 = nodes_list2[[i]]
			seqs2 = nodeIDs_correspondence[which(nodeIDs_correspondence[,1]%in%seqs1),2]
			seqs3 = partition[which(partition[,1]%in%seqs2),2]
			value = as.vector(table(seqs3)[])
			color = names(table(seqs3))
			if (length(value) == 0) value = 1
			if (length(color) == 0) color = "black"
			nodes_values[[i]] = value
			nodes_cols[[i]] = color
		}
	network = graph_from_data_frame(as.data.frame(links[,1:2]), directed=F, vertices=row.names(nodes))
	test.layout = layout_(network, with_dh(weight.edge.lengths=edge_density(network)/scaleFactor2)); par(mar=c(0,0,0,0))
	plot(network, vertex.shape="pie", vertex.pie=nodes_values, vertex.pie.color=nodes_cols, vertex.size=nodes[,"diameter"],
		 layout = test.layout, edge.label=links[,"label"], edge.label.cex=0.5, edge.label.color="black", #edge.label.family="Calibri",
		 edge.label.font=2, vertex.pie.lwd=0.1, vertex.pie.border.lwd=0.1, rescale=T, vertex.label=NA)
	legend("bottomright", colours[,1], text.col="gray30", pch=15, pt.cex=1.5, col=colours[,2], box.lty=0, cex=0.6, y.intersp=1.2)
}

