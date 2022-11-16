
require(Seurat)
require(RColorBrewer)

# background_function: cleanup_data
#
cleanup_data <- function(seurat_slice_obj, query_genes, ctrl_genes=NULL, zonation_score_col="zonation_score", use_sct=TRUE, scale_score=TRUE) {
	if (class(seurat_slice_obj)[1] != "Seurat") {stop("Error: Requires a Seurat Object as input.")}

	# Get expression
	if (use_sct) {
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

# barplot_output == output from zonation_weighted
# annotation_sets == list of vectors where each vector is a set of genes from the same group to be annotated on the plot
# plots lines & labels to annotate groups of genes associated to the same cell-type
add_annotations <- function(barplot_output, annotation_sets) {
	loc <- barplot_output$loc
	hei <- max(barplot_output$height)*0.8
	for (j in 1:length(annotation_sets)) {
		xes <- loc[rownames(loc) %in% annotation_sets[[j]],1]
		if (names(annotation_sets)[j] == "") {next;}
		lines(c(min(xes), max(xes)), c(hei, hei));
		text( x=mean( c(min(xes), max(xes)) ), y=hei, names(annotation_sets)[j], pos=3 )
	}
}


# seurat_slice_obj = Seurat Spatial Transcriptomics object
# query_genes = genes to score for zonation (a vector or list of vectors)
# ctrl_genes = genes to the left of the dashed line used for control purposes (optional)
# zonation_score_col = column of the meta.data to use as the zonation score, or the name of one or more genes to use as the zonation score
# use_sct = use SCT normalization (TRUE) otherwise uses regular @data slot.
# scale_score = whether to rescale the zonation score.
# Description: individual gene expression is rescaled to sum to 100. These are used as weights to calculate the weighted mean zonation score, 
# and weighted standard deviation. Significance is calculated using a linear regression model between the rescaled expression and the zonation score.
#
zonation_weighted <- function(seurat_slice_obj, query_genes, ctrl_genes=NULL, zonation_score_col="zonation_score", use_sct=TRUE, scale_score=TRUE, 
							bar_cols=RColorBrewer::brewer.pal(8, "RdYlBu")) {
	stuff <- cleanup_data(seurat_slice_obj, query_genes, ctrl_genes=NULL, zonation_score_col="zonation_score", use_sct=TRUE, scale_score=TRUE)
	expr_mat <- stuff$expr_mat
	zonation_score <- stuff$zonation_score
	#Error catching
	if (length(zonation_score) == 0) {stop("Zonation score not found.")}

	if (class(query_genes) == "list") {
	# If given a list of gene sets, calculate average expression across all genes in the set.
		calc_zonation_set <- function(gset, expr_mat, zonation_score) {
			# Error catching
			gset <- gset[gset %in% rownames(expr_mat)]
			if (length(gset) < 3) {warning(paste("Warning: only", length(gset), "genes left in the set."))}
	
			if (length(gset) == 1) {
				gexpr <- expr_mat[gset,]
			} else {
				gexpr <- Matrix::colSums(expr_mat[gset,])
			}
			gexpr <- gexpr/sum(gexpr)*100;
			mean <- sum(gexpr*zonation_score)/sum(abs(zonation_score));
			reg <- lm(zonation_score ~ gexpr)
			pval <- summary(reg)$coeff[2,4]
			return(c(mean, pval))
		}
		my_genes <- c(ctrl_genes, query_genes)
		my_names = names(my_genes)
		out <- sapply(my_genes, calc_zonation_set, expr_mat, zonation_score)
		means <- out[1,]
		ps <- out[2,]
	} else {

		# Error catching
		query_genes <- query_genes[query_genes %in% rownames(expr_mat)];
		ctrl_genes <- ctrl_genes[ctrl_genes %in% rownames(expr_mat)];

		if (length(query_genes) == 0) {warning("Warning: No Query Genes found in object")}
		if (length(ctrl_genes) == 0) {warning("No control genes.")}

		# Gather barplot data
		calc_zonation_gene <- function(g, expr_mat, zonation_score) {
			gexpr <- expr_mat[g,]
			gexpr <- gexpr/sum(gexpr)*100;
			mean <- sum(gexpr*zonation_score)/sum(abs(zonation_score));
			reg <- lm(zonation_score ~ gexpr)
			pval <- summary(reg)$coeff[2,4]
			return(c(mean, pval))
		}

		my_genes <- c(ctrl_genes, query_genes)
		my_names = my_genes
		out <- sapply(my_genes, calc_zonation_gene, expr_mat, zonation_score)
		means <- out[1,]
		ps <- out[2,]
	}


	# Plot the result

	# Colour Scale
	data_range <- max(abs(max(means)), abs(min(means)))
	colour_bins <- seq(min(-data_range), max(data_range), length=length(bar_cols)+1)
	colour_bins[1] <- min(1.1*min(means), colour_bins[1]); colour_bins[length(colour_bins)] <- max(1.1*max(means), colour_bins[length(colour_bins)])
	bar_colours <- c(rev(bar_cols)[cut(means, breaks=colour_bins)])

	# Plot Barplot
	locations <- barplot(means, col=bar_colours, ylim=c(min(0,min(means)-0.1*data_range), max(means)+0.1*data_range), 
					ylab="Zonation Score", names=my_names, las=2)
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
	rownames(locations) <- my_names

	# Output
	output <- list(loc=locations, heights=means, pval=ps, ctrl_genes=ctrl_genes, query_genes=query_genes)

	# Lines + Labels
	# If there are multiple labelled sets of genes
	if (! is.null(names(my_genes))  & length(unique(names(my_genes))) != length(my_genes) & length(unique(names(my_genes))) > 1) {
		names(my_genes)[is.na(names(my_genes))] <- ""
		names(my_genes)[my_genes %in% ctrl_genes] <- ""
		groups <- unique(names(my_genes))
		anno_sets <- list()
		for (g in groups) {
			anno_sets[[g]] <- my_genes[names(my_genes) == g]
		}
		add_annotations(output, anno_sets)

	}
	return(output)
}

