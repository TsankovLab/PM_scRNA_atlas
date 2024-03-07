
### Useful functions ###



# Covnert a named list of length of different vectors in a data.frame.
patchvecs = function(vectorList)
	{
	maxLength = max (sapply (vectorList, length))
	tmpL = lapply (vectorList, function(x) c(x, rep (NA, maxLength - length(x))))	
	df = as.data.frame (do.call (cbind,tmpL))
	colnames (df) = names (vectorList)
	return (df)
	}

### Function to create module score of a list gene sets (named), generate a metagroup of highest score across genesets per cell and optionally filter it based on co-expression ###
ModScoreCor = function (seurat_obj, geneset_list, listName, cor_threshold = NULL, pos_threshold = .1, outdir)
        {        
        require (Seurat)	
        message ('Run AddModuleScore')
        seurat_obj = AddModuleScore (seurat_obj, geneset_list)
        seurat_obj@meta.data = seurat_obj@meta.data[, !colnames (seurat_obj@meta.data) %in% names (geneset_list)]
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %in% paste0('Cluster',seq_along(geneset_list))] = names (geneset_list)
        message (paste('Annotate cells based on highest module score and store in column:',paste0(listName, '_r',cor_threshold,'_max')))
        if (length (geneset_list) == 1) 
        	{
        	seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = ifelse (seurat_obj@meta.data[,names (geneset_list)] > pos_threshold, 'pos','neg')
        	pdf (paste0(outdir, listName, '_modulescore_distribution_cor_threshold_',cor_threshold,'_score_',pos_threshold,'.pdf'))
        	hist (seurat_obj@meta.data[,names(geneset_list)])
          abline (v = pos_threshold)
          dev.off()
        	} else {
        	seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = sapply (seq_along(colnames(seurat_obj)), function(x) colnames(seurat_obj@meta.data[,names(geneset_list)])[which.max (seurat_obj@meta.data[x,names(geneset_list)])])        		
        	}
        if (!is.null(cor_threshold))
                {
                message ('cor_threshold provided! Filtering gene sets based on initial correlation to module score')	
                filtered_geneset_list = list()
                geneset_cor_list = list()
                for (i in names(geneset_list))
                        {       
                        geneset_cor = cor (seurat_obj@meta.data[,i], as.matrix(t(seurat_obj@assays$RNA@data[rownames(seurat_obj@assays$RNA@data) %in% geneset_list[[i]],])))
                        geneset_cor_list[[i]] = geneset_cor
                        geneset_cor_names = colnames (geneset_cor)[geneset_cor > cor_threshold]
                        geneset_cor_names = geneset_cor_names[!is.na (geneset_cor_names)]
                        filtered_geneset_list[[i]] = geneset_cor_names
                        }
                if (!is.null (outdir)) 
                        {
                        lapply (seq_along(filtered_geneset_list), function(x) write.csv (filtered_geneset_list[[x]], paste0(outdir,'Corfiltered_Module_score_gene_list_', names(filtered_geneset_list)[x],'.csv')))
                        pdf (paste0(outdir, listName, 'Corfiltered_modulescore_distribution.pdf'))
                        lapply (seq_along(filtered_geneset_list), function(x) 
                                {
                                hist (geneset_cor_list[[x]], title = names(geneset_cor_list)[x])
                                abline (v = cor_threshold)
                                })
                        dev.off()
                        }
                message ('Re-run AddModuleScore using corfiltered genes')
                seurat_obj = AddModuleScore (seurat_obj, filtered_geneset_list)
                seurat_obj@meta.data = seurat_obj@meta.data[, !colnames (seurat_obj@meta.data) %in% paste0(names(geneset_list),'_r',cor_threshold)]
                colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %in% paste0('Cluster',seq_along(geneset_list))] = paste0(names(geneset_list),'_r',cor_threshold)
                if (length (geneset_list) == 1) 
                	{
        					seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = ifelse (seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)] > pos_threshold, 'pos','neg')
        					pdf (paste0(outdir, listName, '_modulescore_distribution_cor_threshold_',cor_threshold,'_score_',pos_threshold,'.pdf'))
        					hist (seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)])
          				abline (v = pos_threshold)
          				dev.off()
        					} else {
                	seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = sapply (seq_along(colnames(seurat_obj)), function(x) colnames(seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)])[which.max (seurat_obj@meta.data[x,paste0(names(geneset_list),'_r',cor_threshold)])])        
                	}
                }
        return (seurat_obj)
        } 



# Compute overlaps of all combinations of vectors within a list
ovmat = function (ovlist, df=FALSE, ov_threshold=0.5, compare_lists= NULL, palette=NULL)
  {
  if (!df) require (ComplexHeatmap)	
  comp_mat = sapply (ovlist, function(x) sapply(ovlist, function (y) sum (unique(x) %in% unique(y)) / min(c(length(unique(x)),length(unique(y))))))
  if (!is.null (compare_lists)) comp_mat = comp_mat[compare_lists[[1]],compare_lists[[2]]]
  if (df) 
  	{
  	return (comp_mat)
  	} else {
    if (is.null(palette)) palette = viridis::mako(100)
  	return (Heatmap (
  		comp_mat, 
  		col=palette,
  		cell_fun = function(j, i, x, y, width, height, fill) {
      if (comp_mat[i,j] > ov_threshold) grid.text (sprintf("%.1f", comp_mat[i, j]), x, y, gp = gpar (fontsize=5))      
      }))      
  	}
  }


# Dotplot showing expression of a gene across two meta groups. Need a Seurat obj
geneDot = function (
  seurat_obj = srt,
	#mat_norm = srt@assays$RNA@data,
	#mat_counts = srt@assays$RNA@counts,
	gene = NULL, 
	x = NULL, # Vector of metagroup1 of equal length to mat columns, If multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
	y = NULL, # vector of metagroup2 of equal length to mat columns
	x_name = 'genes',
  assay='RNA',
	y_name = 'clusters',
	min_expression = 0,
	facet_ncol = 5,
	lim_expression = NULL,
	scale.data = TRUE, # scale data when multiple genes are given
	plotcol = viridis::viridis(100),
	include_NA = TRUE,
  swap_axes = FALSE,
  returnDF = FALSE
	#... # arguments to pass to facet wrap
	)
	{
  require ("scCustomize")
  require ('tidyr')
  if (exists('levels_x')) rm ('levels_x')
  if (exists('levels_y')) rm ('levels_y')
	if (all (grepl ('^\\d', seurat_obj@meta.data[,y]))) seurat_obj@meta.data[,y] = paste0('C',seurat_obj@meta.data[,y])
  if (is.factor (seurat_obj@meta.data[,x])) levels_x = levels (seurat_obj@meta.data[,x]) else
  levels_x = unique (seurat_obj@meta.data[,x])
  if (is.factor (seurat_obj@meta.data[,y])) levels_y = levels (seurat_obj@meta.data[,y]) else 
  levels_y = unique (seurat_obj@meta.data[,y])
  
  seurat_obj@meta.data[,x] = gsub ('[-_ \\+]', '', seurat_obj@meta.data[,x])
  seurat_obj@meta.data[,y] = gsub ('[-_ \\+]', '', seurat_obj@meta.data[,y])
  
  if (exists ('levels_x')) levels_x = gsub ('[-_ \\+]', '', levels_x)
  if (exists ('levels_y')) levels_y = gsub ('[-_ \\+]', '', levels_y)
    
  # Compute percentage expression per cell group
  percent = Percent_Expressing (seurat_object = seurat_obj, assay=assay, threshold = min_expression, features = gene, group_by = x, split_by = y)
  percent$gene = rownames(percent)
  percent = gather (percent, groups, expression, 1:(ncol(percent)-1))
  percent$key = paste0(percent$gene, '_', percent$groups)
  colnames (percent)[colnames(percent) == 'expression'] = 'percent'
  
  # Compute average expression per cell group
  seurat_average = log2 (as.data.frame (AverageExpression (seurat_obj, assay=assay, features= gene, group.by = c(x,y))[[1]]) + 1)
  if (length(gene) == 1) rownames(seurat_average) = gene
  seurat_average$gene = rownames(seurat_average)
  seurat_average = gather (seurat_average, groups, expression, 1:(ncol(seurat_average)-1))
  seurat_average$key = paste0(seurat_average$gene, '_', seurat_average$groups)
  
  # Merge percentage and average expression data.frames
  dot.df = cbind (seurat_average, percent[match (seurat_average$key, percent$key),])
  dot.df = dot.df[,!duplicated (colnames(dot.df))]

  # scale data if scale is TRUE
  if (scale.data) dot.df = transform (dot.df, expression = ave (expression, gene, FUN = function(x) scale(x, scale=T, center=T)))

  dot.df$expression[dot.df$percent == 0] = NA
	dot.df$percent[dot.df$percent == 0] = NA

  dot.df$x_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[1]), levels = levels_x)
  dot.df$y_axis = factor (dot.df$gene, levels = gene) 
  if (!is.null (y)) dot.df$y_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[2]), levels = levels_y)
	if (!is.null (y) & length(gene) > 1) 
    {
    dot.df$y_axis = factor (dot.df$gene, levels = gene) 
    dot.df$z_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[2]), levels = levels_y)
    }

  #dot.df$x_axis = factor (dot.df$x_axis, levels = dot.df$x_axis)
  if (swap_axes)  colnames (dot.df)[match (c('x_axis','y_axis'), colnames (dot.df))] = c('y_axis','x_axis')

	p = ggplot (data = dot.df, aes (x= x_axis, y= y_axis)) +
  	geom_point (shape=21, aes (fill= expression, size = percent), alpha=0.7,colour='black', stroke=0.3) +
    labs (x = x_name, y = y_name, title = ifelse(length(gene) > 1,'',gene), subtitle = paste('Min expression >', min_expression)) +
     scale_shape (solid = FALSE) +
  	theme_classic () +
  	theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1))
  if (!is.null (y) & length(gene) > 1) p = p + facet_wrap (~z_axis)  
  if (length (plotcol) > 3) p = p + scale_fill_gradientn (colors = plotcol) else
    p = p + scale_fill_gradient2 (low = plotcol[1], mid = plotcol[2], high=plotcol[3])	

	if (returnDF) return(dot.df) else 
  return (p)
	}

# Generate barplots or boxplots of meta groups proportions specified
# meta_groups: vector or meta group names:
# 1) first meta_group is the group on which proportions are calculated
# 2) second meta_group split first meta_group on x axes  
#	3) third meta_group will group barplots separately
# if splits include only one value runs barplot instead
cellComp = function (
	seurat_obj = NULL, 
	metaGroups = NULL, # vector of at least 3 metaGroups e.g. c('orig.ident','celltypes','celltypes'),
	plot_as = 'box', # box or bar 
	pal = NULL,
	prop = TRUE,
	ptable_factor = 1, # specify which column of the data.frame or seurat object metadata should be used to compute proportions
	facet_ncol = 20,
	facet_scales = 'free',
	subset_prop = NULL, # subset prop table by any group in any column
	removeNA = TRUE,
	returnDF = FALSE
	)
	{
	require (ggplot2)	
	require (ggpubr)	
	if (is.data.frame (seurat_obj))
		{
		meta_groups_df = seurat_obj[,metaGroups]	
		} else {
		meta_groups_df = seurat_obj@meta.data[,metaGroups]
		}
	# Refactor to remove 0 groups
	#meta_groups_df =  as.data.frame(lapply(unclass(meta_groups_df),as.character),stringsAsFactors=T)
	if(is.null(pal)) pal = rainbow (length(unique(meta_groups_df[,2])))
	#if(is.null(pal) & plot_as == 'box') pal = rainbow (length(unique(meta_groups_df[,3])))
 	if (prop)
 		{
 		ccomp_df = as.data.frame (prop.table (table (meta_groups_df),ptable_factor))
 		ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
 		} else {
 		ccomp_df = as.data.frame (table (meta_groups_df))	
 		}

 	if(removeNA) ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
 	if (!is.null (subset_prop)) 
 		{
 		subset_col = unlist(sapply (seq(ncol(ccomp_df)), function(x) if(any(ccomp_df[,x] %in% subset_prop)) colnames(ccomp_df)[x]))
 		ccomp_df = ccomp_df[ccomp_df[,subset_col] %in% subset_prop,]
 		}
 	#colnames (ccomp_df) = c(paste0('Var_',seq_along(metaGroups)), 'proportion')  
 	if (plot_as == 'box')
 		{
 		p = ggplot (ccomp_df, aes_string (x= metaGroups[2], y= 'Freq')) +    		
  			theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    		scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups > 2)) p = p + geom_boxplot(aes_string (fill= metaGroups[3]),alpha = 0.7, lwd=.2, outlier.shape = NA)
    else p = p + geom_boxplot(aes_string (fill= metaGroups[2]),alpha = 0.7, lwd=.2, outlier.shape = NA)
  	if (length(metaGroups) > 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[4])), scales=facet_scales, ncol=facet_ncol)
  	}		
  if (plot_as == 'bar')
 		{
 		p = ggplot (ccomp_df, aes_string (x= metaGroups[1], y= 'Freq')) +
  			geom_bar(position="stack", stat="identity", aes_string(fill= metaGroups[2])) +
    		#geom_bar(position="dodge", stat="identity", aes_string(fill= metaGroups[2])) +
    		theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    		scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
  	if (length(metaGroups) == 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[3])), scales=facet_scales, ncol=facet_ncol)
  	}
  if (returnDF) return(ccomp_df) else 
  return (p)
	}
	


# Build clustered dotplot with GSEA terms on y axes and clusters on x axes
dotGSEA = function (
	enrichmentsTest_list,
	type = c('fgsea','enrich'),
	top_pathways = NULL, 
	padj_threshold = 0.05,
	cluster_rows = T,
	cluster_cols = T,
	remove_ns_modules = TRUE)
	{
	require (tidyr)	
	if (type == 'fgsea')
		{
		sig_terms = unique(unlist(sapply (enrichmentsTest_list, function(x) x[x$padj < padj_threshold,'pathway'])))
		if (length(sig_terms) == 0) return (NULL)
		if (!is.null(top_pathways)) sig_terms = unique(unlist(sapply (enrichmentsTest_list, function(x) 
				{
				x =	na.omit(x)
				x = x[x$padj < padj_threshold,]
				x_top = x[order(-x$NES),]
				x_top = head(x_top[order(x_top$padj),'pathway'], top_pathways)
				x_bot = x[order(x$NES),]
				x_bot = head(x_bot[order(-x_bot$padj),'pathway'], top_pathways)
				c(x_top,x_bot)
				})))
		if (length(sig_terms) == 1)
			{
			mat_pvalue = t(as.matrix (sapply (enrichmentsTest_list, function(x) x$pval[match(sig_terms, x$pathway)])))
			mat_sizeLog = t(as.matrix (sapply (enrichmentsTest_list, function(x) x$NES[match(sig_terms, x$pathway)])))	
			} else {
			mat_pvalue = sapply (enrichmentsTest_list, function(x) x$pval[match(sig_terms, x$pathway)])
			mat_sizeLog = sapply (enrichmentsTest_list, function(x) x$NES[match(sig_terms, x$pathway)])	
			}
		}
	if (type == 'enrich')
		{
		sig_terms = na.omit(unique(unlist(sapply (enrichmentsTest_list, function(x) x[x$p.adjust < padj_threshold,'ID']))))
		if (length(sig_terms) == 0) return (NULL)
		if (!is.null(top_pathways)) sig_terms = unique(as.vector(unlist (sapply (enrichmentsTest_list, function(x) na.omit(head(x[x$p.adjust < padj_threshold,'ID'], top_pathways))))))

		if (length(sig_terms) == 1)
			{
			mat_pvalue = t(as.matrix(sapply (enrichmentsTest_list, function(x) x$p.adjust[match(sig_terms, x$ID)])))
			mat_sizeLog = t(as.matrix(sapply (enrichmentsTest_list, function(x) x$Count[match(sig_terms, x$ID)])))	
			} else {
			mat_pvalue = sapply (enrichmentsTest_list, function(x) x$p.adjust[match(sig_terms, x$ID)])
			mat_sizeLog = sapply (enrichmentsTest_list, function(x) x$Count[match(sig_terms, x$ID)])	
			}
		}
	
	mat_pvalue[is.na(mat_pvalue)] = 1
	mat_sizeLog[is.na(mat_sizeLog)] = 0
	
rownames (mat_pvalue) = sig_terms
colnames (mat_pvalue) = names (enrichmentsTest_list)
rownames (mat_sizeLog) = sig_terms
colnames (mat_sizeLog) = names (enrichmentsTest_list)
mat_sizeLog = as.data.frame (mat_sizeLog)
mat_pvalue = as.data.frame (mat_pvalue)

	if(cluster_rows & min (dim(mat_sizeLog)[1]) >= 2) 
		{
		d_row = dist (mat_sizeLog)
		d_row[is.na(d_row)] = max (d_row)	
		hc1_row = hclust (d_row, method = "ward.D")
		mat_sizeLog = mat_sizeLog[hc1_row$order,]
		mat_pvalue = mat_pvalue[hc1_row$order,]
		}
	if(cluster_cols & min (dim(mat_sizeLog)[2]) >= 2) 
		{
		d_col = dist (t(mat_sizeLog))	
		d_col[is.na(d_col)] = max (d_col)
		hc1_col = hclust (d_col, method = "ward.D")
		mat_sizeLog = mat_sizeLog[,hc1_col$order]
		mat_pvalue = mat_pvalue[,hc1_col$order]
		}
#mat_sizeLog = as.data.frame (mat_sizeLog)
mat_sizeLog[mat_pvalue > padj_threshold] = NA			
mat_pvalue[mat_pvalue > padj_threshold] = NA
mat_pvalue = -log10(as.data.frame (mat_pvalue))	
#colnames (mat_pvalue) = names (enrichmentsTest_list)

mat_pvalue$pathway = rownames (mat_pvalue)
if (remove_ns_modules) # remove modules with no enriched pathways
	{
	remove_cols	= apply (mat_pvalue, 2, function(x) !all (is.na(x)))
	mat_pvalue = mat_pvalue[, remove_cols, drop=F]
	mat_sizeLog = mat_sizeLog[, colnames(mat_pvalue)[colnames(mat_pvalue)!='pathway'], drop=F]
	}
if (type == 'fgsea')
	{
	gsea_df_pvalue = gather (mat_pvalue, cluster, nlogPvalue, colnames(mat_pvalue)[1]:colnames(mat_pvalue)[ncol(mat_pvalue)-1], factor_key=TRUE)
	gsea_df_NES = gather (mat_sizeLog, cluster2, NES, colnames(mat_sizeLog)[1]:colnames(mat_sizeLog)[ncol(mat_sizeLog)], factor_key=TRUE)
	gsea_df = cbind (gsea_df_NES, gsea_df_pvalue)
	#gsea_df$cluster = factor (gsea_df$cluster, levels = unique(gsea_df$cluster))
	gsea_df$pathway = factor (gsea_df$pathway, levels = unique(gsea_df$pathway))
	p = ggplot (data = gsea_df, aes (x=cluster, y=pathway,
	      fill= NES, size=nlogPvalue)) +
	      geom_point (shape=21, color='black') +
	      #scale_fill_gradient2 (low = 'blue',mid='white', high = 'red', midpoint=0) +
	      scale_fill_gradient2(
      	low = "#67A9CF", 
      	mid = "#F7F7F7", 
      	high = "#EF8A62", 
      	midpoint = 0
    		) +
	      labs (x = 'cluster', y = '-log10(pvalue)') +
	      theme_minimal() + 
	      theme(text = element_text(size=11)) +
	      theme(axis.text.y = element_text (angle = 0, vjust = 0.5, hjust=1)) +
	      theme(axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1))
	return (p)	
	}
if (type == 'enrich')
	{	
	enrich_df_pvalue = gather (mat_pvalue, cluster, nlogPvalue, colnames(mat_pvalue)[1]:colnames(mat_pvalue)[ncol(mat_pvalue)-1], factor_key=TRUE)
	enrich_df_Count = gather (mat_sizeLog, cluster2, Count, colnames(mat_sizeLog)[1]:colnames(mat_sizeLog)[ncol(mat_sizeLog)], factor_key=TRUE)
	enrich_df = cbind (enrich_df_Count, enrich_df_pvalue)
	#gsea_df$cluster = factor (gsea_df$cluster, levels = unique(gsea_df$cluster))
	enrich_df$pathway = factor (enrich_df$pathway, levels = unique(enrich_df$pathway))
	p = ggplot (data = enrich_df, aes (x=cluster, y=pathway,
	      fill= nlogPvalue, size=Count)) +
	      geom_point (shape=21, color='black') +
	      scale_fill_distiller (palette = "Spectral") +
	      labs (x = 'cluster', y = '-log10(pvalue)') +
	      theme_minimal() + 
	      theme(text = element_text(size=11)) +
	      theme(axis.text.y = element_text (angle = 0, vjust = 0.5, hjust=1)) +
	      theme(axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1))
	return (p)		
	}
	}

