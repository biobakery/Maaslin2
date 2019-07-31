###############################################################################
# MaAsLin2 visualizations

# Copyright (c) 2018 Harvard School of Public Health

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################

# Author: Gholamali (Ali) rahnavard
# Email: gholamali.rahnavard@gmail.com
# This script includes functions for visualizing overall output of MaAsLin2 and
# individual associations as scatterplot and boxplot

# Load libararies
for (lib in c('ggplot2', "grid", 'pheatmap')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# MaAsLin2 theme based on Nature journal requirements
nature_theme <- function(x_axis_labels, y_label) {
    # set default text format based on categorical and length
    angle = NULL
    hjust = NULL
    size = 8
    if (max(nchar(x_axis_labels), na.rm=TRUE) > 5) {
        angle = 45
        hjust = 1
        size = 6
    }
    axis_title_size = 10
    if (nchar(y_label) > 15) {
        axis_title_size = 8
    }
    if (nchar(y_label) > 25) {
        axis_title_size = 6
    }
    return ( ggplot2::theme_bw() + ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = size, vjust = 1, hjust = hjust, angle = angle),
        axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
        axis.title = ggplot2::element_text(size = axis_title_size),
        plot.title = ggplot2::element_text(size = 7, face = 'bold'),
        legend.title = ggplot2::element_text(size = 6, face = 'bold'),
        legend.text = ggplot2::element_text(size = 6),
        axis.line = ggplot2::element_line(colour = 'black', size = .25),
        axis.line.x = ggplot2::element_line(colour = 'black', size = .25),
        axis.line.y = ggplot2::element_line(colour = 'black', size = .25),
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank())
   )
}


# MaAsLin2 heatmap function for overall view of associations
maaslin2_heatmap <-
    function(
        output_results,
        title = NA,
        cell_value = 'qval',
        data_label = 'data',
        metadata_label = 'metadata',
        border_color = 'grey93',
        color = colorRampPalette(c("darkblue", "grey90", "darkred")),
        col_rotate = 90,
        first_n = 50) {

        # read MaAsLin output
        df <- read.table(
            output_results,
            header = TRUE,
            sep = "\t",
            fill = TRUE,
            comment.char = "" ,
            check.names = FALSE
        )

        title_additional <- ""
        
        title_additional <- ""
        if (!is.na(first_n) & first_n > 0 & first_n < dim(df)[1]) {
            if (cell_value == 'coef') {
                df <- df[order(-abs(df[cell_value])) , ]
            } else{
                df <- df[order(df[cell_value]), ]
            }
            # get the top n features with significant associations
            df_sub <- df[1:first_n,]
            for (first_n_index in seq(first_n, dim(df)[1]))
            {
                if (length(unique(df_sub$feature)) == first_n)
                {
                    break
                }
                df_sub <- df[1:first_n_index,]
            }
            # get all rows that have the top N features
            df <- df[which(df$feature %in% df_sub$feature),]
            title_additional <- paste("Top", first_n, sep=" ")
        }
        
        if (dim(df)[1] < 2) {
            print('There are no associations to plot!')
            return(NULL)
        }
        
        metadata <- df$metadata
        data <- df$feature
        value <- NA
        
        # values to use for coloring the heatmap
        # and set the colorbar boundaries
        if (cell_value == "pval") {
            value <- -log(df$pval) * sign(df$coef)
            value <- pmax(-20, pmin(20, value))
            if (is.null(title))
                title <- "(-log(pval)*sign(coeff))"
        } else if (cell_value == "qval") {
            value <- -log(df$qval) * sign(df$coef)
            value <- pmax(-20, pmin(20, value))
            if (is.null(title))
                title <- "(-log(qval)*sign(coeff))"
        } else if (cell_value == "coef") {
            value <- df$coef
            if (is.null(title))
                title <- "(coeff)"
        }
   
        if (title_additional!="") {
            title <- paste(title_additional, "features with significant associations", title, sep=" ")
        } else {
            title <- paste("Significant associations", title, sep=" ")
        }

        n <- length(unique(data))
        m <- length(unique(metadata))
        
        if (n < 2) {
            print(
                paste(
                    "There is not enough features in the associations",
                    "to create a heatmap plot.",
                    "Please review the associations in text output file.")
            )
            return(NULL)
        }
        
        if (m < 2) {
            print(
                paste(
                    "There is not enough metadata in the associations",
                    "to create a heatmap plot.",
                    "Please review the associations in text output file.")
            )
            return(NULL)
        }
        
        a = matrix(0, nrow = n, ncol = m)
        a <- as.data.frame(a)
        
        rownames(a) <- unique(data)
        colnames(a) <- unique(metadata)
        for (i in seq_len(dim(df)[1])) {
            if (abs(a[as.character(data[i]), 
                    as.character(metadata[i])]) > abs(value[i]))
                next
            a[as.character(data[i]), as.character(metadata[i])] <- value[i]
        }
       
        # get the range for the colorbar
        max_value <- ceiling(max(a))
        min_value <- ceiling(min(a))
        range_value <- max(c(abs(max_value),abs(min_value)))
        breaks <- seq(-1*range_value, range_value, by = 1)

        p <- NULL
        tryCatch({
            p <-
                pheatmap::pheatmap(
                    a,
                    cellwidth = 5,
                    cellheight = 5,
                    # changed to 3
                    main = title,
                    fontsize = 6,
                    kmeans_k = NA,
                    border = TRUE,
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    scale = "none",
                    cluster_rows = FALSE,
                    cluster_cols = TRUE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    legend = TRUE,
                    border_color = border_color,
                    color = color(range_value*2),
                    breaks = breaks,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    display_numbers = matrix(ifelse(
                        a > 0.0, "+", ifelse(a < 0.0, "-", "")), nrow(a))
                )
        }, error = function(err) {
            logging::logerror("Unable to plot heatmap")
            logging::logerror(err)
        })
        return(p)
    }

save_heatmap <-
    function(
        results_file,
        heatmap_file,
        title = NULL,
        cell_value = "qval",
        data_label = 'data',
        metadata_label = 'metadata',
        border_color = "grey93",
        color = colorRampPalette(c("blue", "grey90", "red")),
        first_n = 50) {

        # generate a heatmap and save it to a pdf
        heatmap <-
            maaslin2_heatmap(
                results_file,
                title,
                cell_value,
                data_label,
                metadata_label,
                border_color,
                color,
                first_n)
        
        if (!is.null(heatmap)) {
            pdf(heatmap_file)
            print(heatmap)
            dev.off()
        }
        
    }

maaslin2_association_plots <-
    function(
        metadata,
        features,
        output_results,
        write_to = './')
    {
        #MaAslin2 scatter plot function and theme
        
        # combine the data and metadata to one datframe using common rows
        # read MaAsLin output
        if (is.character(metadata)) {
            metadata <- read.table(
                metadata,
                header = TRUE,
                row.names = 1,
                sep = "\t",
                fill = FALSE,
                comment.char = "" ,
                check.names = FALSE
            )
        }
        if (is.character(features)) {
            features <- read.table(
                features,
                header = TRUE,
                row.names = 1,
                sep = "\t",
                fill = FALSE,
                comment.char = "" ,
                check.names = FALSE
            )
        }
        
        common_rows <- intersect(rownames(features), rownames(metadata))
        input_df_all <-
            cbind(features[common_rows, , drop = FALSE], 
                metadata[common_rows, , drop = FALSE])
        
        # read MaAsLin output
        if (is.character(output_results)) {
            output_df_all <- read.table(
                output_results,
                header = TRUE,
                row.names = NULL,
                sep = "\t",
                fill = FALSE,
                comment.char = "" ,
                check.names = FALSE
            )
        } else {
            output_df_all <- output_results
        }
        
        if (dim(output_df_all)[1] < 1) {
            print('There are no associations to plot!')
            return(NULL)
        }
        
        logging::loginfo(
            paste("Plotting associations from most",
                "to least significant,",
                "grouped by metadata"))
        metadata_types <- unlist(output_df_all[, 'metadata'])
        metadata_labels <-
            unlist(metadata_types[!duplicated(metadata_types)])
        metadata_number <- 1
        
        for (label in metadata_labels) {
            # for file name replace any non alphanumeric with underscore
            plot_file <-
                paste(
                    write_to,
                    "/",
                    gsub("[^[:alnum:]_]", "_", label),
                    ".pdf",
                    sep = "")
            data_index <- which(label == metadata_types)
            logging::loginfo("Plotting data for metadata number %s, %s",
                metadata_number,
                label)
            pdf(
                plot_file,
                width = 2.65,
                height = 2.5,
                onefile = TRUE)

            x <- NULL
            y <- NULL 
            for (i in data_index) {
                x_label <- as.character(output_df_all[i, 'metadata'])
                y_label <- as.character(output_df_all[i, 'feature'])
                qval <- as.numeric(output_df_all[i, 'qval'])
                coef_val <- as.numeric(output_df_all[i, 'coef'])
                input_df <- input_df_all[c(x_label, y_label)]
                colnames(input_df) <- c("x", "y")
                
                # if Metadata is continuous generate a scatter plot
                # Continuous is defined as numerical with more than 
                # 2 values (to exclude binary data)
                temp_plot <- NULL
                if (is.numeric(input_df[1, 'x']) &
                        length(unique(input_df[['x']])) > 2) {
                    logging::loginfo(
                        "Creating scatter plot for continuous data, %s vs %s",
                        x_label,
                        y_label)
                    temp_plot <- ggplot2::ggplot(
                        data = input_df,
                            ggplot2::aes(
                                as.numeric(as.character(x)),
                                as.numeric(as.character(y))
                                )) +
                            ggplot2::geom_point(
                                fill = 'darkolivegreen4',
                                color = 'black',
                                alpha = .5,
                                shape = 21,
                                size = 1,
                                stroke = 0.15
                            ) +
                            ggplot2::scale_x_continuous(
                                limits = c(min(
                                    input_df['x']), max(input_df['x']))) +
                            ggplot2::scale_y_continuous(
                                limits = c(min(
                                    input_df['y']), max(input_df['y']))) +
                            ggplot2::stat_smooth(
                                method = "glm",
                                size = 0.5,
                                color = 'blue',
                                na.rm = TRUE
                            ) +
                            ggplot2::guides(alpha = 'none') + 
                            ggplot2::labs("") +
                            ggplot2::xlab(x_label) + 
                            ggplot2::ylab(y_label) + 
                            nature_theme(input_df[, 'x'], y_label) +
                            ggplot2::annotate(
                                geom = "text",
                                x = Inf,
                                y = Inf,
                                hjust = 1,
                                vjust = 1,
                                label = sprintf(
                                    "FDR: %s\nCoefficient: %s\nN: %s",
                                    formatC(qval, format = "e", digits = 3),
                                    formatC(coef_val, format = "e", digits = 2),
                                    formatC(length(input_df[, 'x']))
                                ) ,
                                color = "black",
                                size = 2,
                                fontface = "italic"
                            )
                } else{
                    # if Metadata is categorical generate a boxplot
                    ### check if the variable is categorical
                    
                    logging::loginfo(
                        "Creating boxplot for catgorical data, %s vs %s",
                        x_label,
                        y_label)
                    input_df['x'] <- lapply(input_df['x'], as.character)

                    # count the Ns for each group
                    x_axis_label_names <- unique(input_df[['x']])
                    for (name in x_axis_label_names) {
                        total <- length(which(input_df[['x']] == name))
                        new_n <- paste(name, " (n=", total, ")", sep="")
                        input_df[which(input_df[['x']] == name),'x'] <- new_n
                    }
                    temp_plot <-
                        ggplot2::ggplot(
                            data = input_df, ggplot2::aes(factor(x), y)) +
                        ggplot2::geom_boxplot(
                            ggplot2::aes(fill = x),
                            outlier.alpha = 0.0,
                            na.rm = TRUE,
                            alpha = .5,
                            show.legend = FALSE
                        ) +
                        ggplot2::geom_point(
                            ggplot2::aes(fill = x),
                            alpha = 0.75 ,
                            size = 1,
                            shape = 21,
                            stroke = 0.15,
                            color = 'black',
                            position = ggplot2::position_jitterdodge()
                        ) +
                        ggplot2::scale_fill_brewer(palette = "Spectral")
                    
                    # format the figure to default nature format
                    # remove legend, add x/y labels
                    temp_plot <- temp_plot + 
                        nature_theme(input_df[, 'x'], y_label) +
                        ggplot2::theme(
                            panel.grid.major = ggplot2::element_blank(),
                            panel.grid.minor = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank(),
                            axis.line = ggplot2::element_line(colour = "black")
                        ) +
                        ggplot2::xlab(x_label) +
                        ggplot2::ylab(y_label) +
                        ggplot2::theme(legend.position = "none") +
                        ggplot2::annotate(
                            geom = "text",
                            x = Inf,
                            y = Inf,
                            hjust = 1,
                            vjust = 1,
                            label = sprintf(
                                "FDR: %s\nCoefficient: %s",
                                formatC(qval, format = "e", digits = 3),
                                formatC(coef_val, format = "e", digits = 2)
                            ) ,
                            color = "black",
                            size = 2,
                            fontface = "italic"
                        )
                }
                stdout <- capture.output(print(temp_plot), type = "message")
                if (length(stdout) > 0)
                    logging::logdebug(stdout)
            }
            dev.off()
            metadata_number <- metadata_number + 1
        }
    }
