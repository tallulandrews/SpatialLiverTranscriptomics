
require(Seurat)
require(RColorBrewer)

# background_function: cleanup_data
#
cleanup_data <- function(seurat_slice_obj, query_genes, ctrl_genes=NULL, zonation_score_col="zonation_score", use_sct=TRUE, scale_score=TRUE) {
	if (class(seurat_slice_obj)[1] != "Seurat") {stop("Error: Requires a Seurat Object as input.")}

	# Get expression
	if (use.sct) {
		expr_mat <- GetAssayData(GetAssay(seurat_slice_obj, "SCT"), slot="data")
	} else {
		expr_mat <- GetAssayData(GetAssay(seurat_slice_obj, "Spatial"), slot="data")
	}


	# Zonation Score
	zonation_score <- c();
	if (zonation_score_col %in% colnames(seurat_slice_obj@meta.data)) {
		zonation_score <- seurat_slice_obj@meta.data[,zonation_score_col]
	}
	if (zonation_score_col %in% rownames(seurat_slice_obj)) {
		if (sum(zonation_score_col %in% rownames(seurat_slice_obj)) == 1) {
			zonation_score <- expr_mat[zonation_score_col,]
		} else {
			zonation_score <- Matrix::colSums(expr_mat[zonation_score_col,])
		}
	}
	if (scale_score) {
		zonation_score <- scale(zonation_score)
	}
	return(list(expr_mat=expr_mat, zonation_score=zonation_score))
}



# seurat_slice_obj = Seurat Spatial Transcriptomics object
# query_genes = genes to score for zonation
# ctrl_genes = genes to the left of the dashed line used for control purposes (optional)
# zonation_score_col = column of the meta.data to use as the zonation score, or the name of one or more genes to use as the zonation score
# use_sct = use SCT normalization (TRUE) otherwise uses regular @data slot.
# scale_score = whether to rescale the zonation score.
# Description: individual gene expression is rescaled to sum to 100. These are used as weights to calculate the weighted mean zonation score, 
# and weighted standard deviation. Significance is calculated using a linear regression model between the rescaled expression and the zonation score.
#
zonation_genes_weighted <- function(seurat_slice_obj, query_genes, ctrl_genes=NULL, zonation_score_col="zonation_score", use_sct=TRUE, scale_score=TRUE) {
	violin_colours <- RColorBrewer::brewer.pal(8, "RdYlBu")
	stuff <- cleanup_data(seurat_slice_obj, query_genes, ctrl_genes=NULL, zonation_score_col="zonation_score", use_sct=TRUE, scale_score=TRUE)
	expr_mat <- stuff$expr_mat
	zonation_score <- stuff$zonation_score

	# Error catching
	query_genes <- query_genes[query_genes %in% rownames(expr_mat)];
	ctrl_genes <- ctrl_genes[ctrl_genes %in% rownames(expr_mat)];

	if (length(query_genes) == 0) {warning("Warning: No Query Genes found in object")}
	if (length(ctrl_genes) == 0) {warning("No control genes.")}
	if (length(zonation_score) == 0) {stop("Zonation score not found.")}

	# Gather barplot data

	means <- c()
	ps <- c()
	
	my_genes <- c(ctrl_genes, query_genes)
	for (g in my_genes) {
		expr <- expr_mat[g,]
		expr <- expr/sum(expr)*100;
		means <- c(means, mean(expr*zonation_score));
		reg <- lm(zonation_score ~ expr)
		ps <- c(ps, summary(reg)$coeff[2,4])
	}

	# Plot the result

	# Colour Scale
	data_range <- max(abs(max(means)), abs(min(means)))
	colour_bins <- seq(min(-data_range), max(data_range), length=length(violin_colours)+1)
	colour_bins[1] <- min(1.1*min(means), colour_bins[1]); colour_bins[length(colour_bins)] <- max(1.1*max(means), colour_bins[length(colour_bins)])
	bar_colours <- c(rev(violin_colours)[cut(means, breaks=colour_bins)])

	# Plot Barplot
	locations <- barplot(means, col=bar_colours, ylim=c(min(0,min(means)-0.1*data_range), max(means)+0.1*data_range), 
					ylab="Zonation Score", names=my_genes, las=2)
	label_pos <- rep(3, length(means)); label_pos[means > 0] <- 1;

	# Ctrl line
	if (length(ctrl_genes) > 0) {
		abline(v=mean(locations[c(length(ctrl_genes), length(ctrl_genes)+1)]), lty=3)
	}

	# Significance
	p_text <- signif(ps[1:length(ps)], digits=1)
	p_stars <- rep("", length(p_text));
	p_stars[p_text < 0.01] <- "*"
	p_stars[p_text < 10^-10] <- "**"
	p_stars[p_text < 10^-100] <- "***"
	text(locations, means, p_stars, pos=abs(4-label_pos))
	
	# Output
	rownames(locations) <- my_genes
	return(list(loc=locations, heights=means, pval=ps, ctrl_genes=ctrl_genes, query_genes=query_genes))
}


# seurat_slice_obj = Seurat Spatial Transcriptomics object
# gene_lists = a list of vectors of marker genes for each cell-type, these must be named for what cell-type they represent
# ctrl_types = vector of cell-types to consider controls.
# zonation_score_col = column of the meta.data to use as the zonation score, or the name of one or more genes to use as the zonation score
# use_sct = use SCT normalization (TRUE) otherwise uses regular @data slot.
# scale_score = whether to rescale the zonation score.
# Description: individual gene expression is rescaled to sum to 100. These are used as weights to calculate the weighted mean zonation score, 
# and weighted standard deviation. Significance is calculated using a linear regression model between the rescaled expression and the zonation score.
#
aggregate_score <- function(seurat_slice_obj, gene_lists, ctrl_types=NULL, zonation_score_col="zonation_score", use_sct=TRUE, scale_score=TRUE) {
	violin_colours <- RColorBrewer::brewer.pal(8, "RdYlBu")
	stuff <- cleanup_data(seurat_slice_obj, query_genes, ctrl_genes=NULL, zonation_score_col="zonation_score", use_sct=TRUE, scale_score=TRUE)
	expr_mat <- stuff$expr_mat
	zonation_score <- stuff$zonation_score
	if (length(zonation_score) == 0) {stop("Zonation score not found.")}

	# Error catching
	gene_lists <- lapply(gene_lists, function(x){x[x %in% rownames(seurat_slice_obj)]})
	exclude <- sapply(gene_lists, function(x){length(x)==0})
	if (sum(exclude) > 0) {	
		gene_lists <- gene_lists[-1*which(exclude)]
	}
	ctrl_types <- ctrl_types[ctrl_types %in% names(gene_lists)]

	# Gather barplot data

	means <- c()
	ps <- c()
	
	my_types <- names(my_gene_lists)
	if (!is.null(ctrl_types)) {
		my_types <- c(ctrl_types, my_types[ !(my_types %in% ctrl_types) ])
	}
	for (i in my_types) {
		g <- my_gene_lists[[i]]
		g <- g[g %in% rownames(expr_mat)]
		if (length(g) == 1) {
			expr <- expr_mat[g,]
		} else {
			expr <- Matrix::colSums(expr_mat[g,])
		}
		expr <- expr/sum(expr)*100;
		means <- c(means, mean(expr*zonation_score));
		reg <- lm(zonation_score ~ expr)
		ps <- c(ps, summary(reg)$coeff[2,4])
	}

	# Plot the result

	# Colour Scale
	data_range <- max(abs(max(means)), abs(min(means)))
	colour_bins <- seq(min(-data_range), max(data_range), length=length(violin_colours)+1)
	colour_bins[1] <- min(1.1*min(means), colour_bins[1]); colour_bins[length(colour_bins)] <- max(1.1*max(means), colour_bins[length(colour_bins)])
	bar_colours <- c(rev(violin_colours)[cut(means, breaks=colour_bins)])

	# Plot Barplot
	locations <- barplot(means, col=bar_colours, ylim=c(min(0,min(means)-0.1*data_range), max(means)+0.1*data_range), 
					ylab="Zonation Score", names=my_types, las=2)
	label_pos <- rep(3, length(means)); label_pos[means > 0] <- 1;

	# Ctrl line
	if (length(ctrl_types) > 0) {
		abline(v=mean(locations[c(length(ctrl_types), length(ctrl_types)+1)]), lty=3)
	}

	# Significance
	p_text <- signif(ps[1:length(ps)], digits=1)
	p_stars <- rep("", length(p_text));
	p_stars[p_text < 0.01] <- "*"
	p_stars[p_text < 10^-10] <- "**"
	p_stars[p_text < 10^-100] <- "***"
	text(locations, means, p_stars, pos=abs(4-label_pos))
	
	# Output
	rownames(locations) <- my_types
	return(list(loc=locations, heights=means, pval=ps, ctrl_genes=ctrl_types, gene_lists=gene_lists))
}

add_annotations <- function(barplot_output, annotation_sets) {
	loc <- barplot_output$loc
	hei <- max(barplot_output$height)*0.8
	for (j in 1:length(annotation_sets)) {
		xes <- loc[rownames(loc) %in% annotation_sets[[j]],1]
		lines(c(min(xes), max(xes)), c(hei, hei));
		text(x=mean(c(min(xes), max(xes)), y=hei, names(annotation_sets)[j], pos=3)
	}
}
