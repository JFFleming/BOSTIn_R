library(shiny)
library(phangorn)
library(Biostrings)
library(reticulate)


# Define server logic to read selected file
server <- function(input, output) {
####THIS IS TO DEAL WITH SITE SAT UI INPUT THAT VARIES BASED ON DATA TYPE #####
  output$dynamic_tabs <- renderUI({
    if (input$d_type == "DNA") {
      tabsetPanel(
        tabPanel("C-Score", tableOutput("freq_table")),
        tabPanel("taxon C-Score", tableOutput("taxa_table")),
        tabPanel("taxon C-Score Histogram", plotOutput("sathist")),
        tabPanel("Bostin", uiOutput("sitesat_bostin"))
      )
    } else {
      tabsetPanel(
        tabPanel("DE-Score", tableOutput("freq_table")),
        tabPanel("tDE-Score", tableOutput("taxa_table")),
        tabPanel("tDE-Score Histogram", plotOutput("sathist")),
        tabPanel("Bostin", uiOutput("sitesat_bostin"))
      )
    }
  })
####THIS IS THE BRANCH HETEROGENEITY PORTION####
  observeEvent(input$do,{
    req(input$file1)  # Make sure file is uploaded
    req(input$branchhet)
    
    tryCatch({
      input_fasta <- read.phyDat(input$file1$datapath, format = "fasta", type = input$d_type)
      file_data_blh <<- calculate_LBi(input_fasta)
      
    }, error = function(e) {
      stop(safeError(e))
    })
  
####BLH Output Tabs###
  # Render plot (for example, plot the phylogenetic tree)
  output$treeplot <- renderPlot({
    req(file_data_blh)  # Ensure data is available
    plot(file_data_blh$treeNJ, main = "Phylogenetic Tree")
  })
  #Render Histogram
  output$histplot <- renderPlot({
  req(file_data_blh)  # Ensure data is available
  result_matrix <- file_data_blh$result_matrix  # Access the result matrix
  
  # Extract the upper quartile values
#  UQ <- result_matrix[,2][result_matrix[,2] >= quantile(result_matrix[,2])[4]]
  LB_All <- result_matrix[,2]
  # Plot histogram of upper quartile scores
  hist(LB_All,
       main = "Histogram of LB-Scores",
       xlab = "LB-Scores",
       col = "lightblue",
       border = "black")  # Adjust breaks as necessary
})
  # Render summary
output$summary <- renderTable({
  req(file_data_blh)  # Ensure data is available
  result_matrix <- file_data_blh$result_matrix  # Access the result matrix
  
  # Calculate statistics
  std_dev <- sd(result_matrix[, 2], na.rm = TRUE)
  mean_val <- mean(result_matrix[, 2], na.rm = TRUE)
  min_val <- quantile(result_matrix[, 2], 0)[1]
  lower_quartile <- quantile(result_matrix[, 2], 0.25)[1]
  median_val <- quantile(result_matrix[, 2], 0.5)[1]
  upper_quartile <- quantile(result_matrix[, 2], 0.75)[1]
  max_val <- quantile(result_matrix[, 2], 1)[1]
  
  # Create a data frame for summary statistics
  summary_df <- data.frame(
    Statistic = c("Standard Deviation", "Mean", "Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum"),
    Value = c(std_dev, mean_val, min_val, lower_quartile, median_val, upper_quartile, max_val)
  )
  
  # Return the data frame directly as a table
  return(summary_df)
})
  # Render table for the LongBranch tab
  output$list <- renderTable({
    req(file_data_blh)  # Ensure data is available
    file_data_blh$result_matrix
  }, rownames=TRUE)
  

output$lb_bostin <- renderUI({
  req(file_data_blh)  # Ensure data is available
  result_matrix <- file_data_blh$result_matrix  # Access the result matrix
  
  # Extract statistics from result_matrix
  mean_val <- mean(result_matrix[, 2], na.rm = TRUE)
  std_dev <- sd(result_matrix[, 2], na.rm = TRUE)
  min_val <- quantile(result_matrix[, 2], 0)[1]
  max_val <- quantile(result_matrix[, 2], 1)[1]
  upper_quartile <- quantile(result_matrix[, 2], 0.75)[1]
  lower_quartile <- quantile(result_matrix[, 2], 0.25)[1]
  UQStdDev <- sd(result_matrix[, 2], na.rm = TRUE)  # Assuming stdDevUpperQuartile is the same as std_dev
  LB_MAD <- mad(result_matrix[, 2])
  LB_MED <- median(result_matrix[, 2])
  
  # Calculate the upper quartile boundaries
  UQBound1 <- LB_MED + (2 * LB_MAD)
  UQBound2 <- LB_MED + (3 * LB_MAD)
#  print(UQBound1)
#  print(UQBound2)
  # Flag taxa based on their values
  yellowFlag <- vector()
  redFlag <- vector()
  
  # Loop through taxa and flag them
  for (i in 1:nrow(result_matrix)) {
    if (result_matrix[i, 2] > UQBound2) {
      redFlag <- c(redFlag, rownames(result_matrix)[i])
    }
    if (result_matrix[i, 2] > UQBound1) {
      yellowFlag <- c(yellowFlag, rownames(result_matrix)[i])
    }
  }
  # Construct the old perl output text
output_text <- paste0(
    
    "<p><strong>LB Score</strong></p>",
    
    "<p>To identify long branched taxa, BOSTIn uses the LB-score. The LB-Score measures the percentage deviation of each taxon from the average patristic distance, and so is independent of the actual topology of the tree itself, making it quite useful to identify long branches. BostIn rapidly generates a Neighbour-Joining tree to calculate the LB-Score. This produces an LB-Score that is normally significantly similar, even under large amounts of Long Branch Attraction, but it won't be as accurate as an LB-Score generated under the best possible model. For the purposes of defining the sextile of taxa most likely to cause a long branch attraction artifact, however, it ought to suffice.</p>",
    
    "<p>To read more about this, see the BOSTIn manuscript when it appears in pre-print (I'll add a reference here later!).</p>",
    
    "<p>The taxa specific LB-Scores in your dataset range from ", min_val, " to ", max_val, 
    ", with a median of ", LB_MED, " a mean of ", mean_val, " and a standard deviation of ", std_dev, ".</p>",
    
    "<p>Using the median and the median absolute deviation (similar to a standard deviation, but for medians), we can use the LB-Score to more robustly identify suspect long-branched taxa by assessing which taxa are outside of two, and then 3 median absolute deviations of the median. This is because it is the extremes of branch length heterogeneity that can cause the greatest problems. The median absolute deviation is ", LB_MAD, ", and so the median plus two median absolute deviations is ", UQBound1, ", while plus three median absolute deviations ", UQBound2, ".</p>",
    
    "<p>We've identified taxa beyond these bounds as yellow flags and red flags respectively, as with the other measurements in BOSTIn.</p>",
    
    "<p><strong>Your Red Flag Taxa are:</strong></p>",
    "<ul>", paste("<li>", redFlag, "</li>", collapse = ""), "</ul>",

    "<p><strong>Your Yellow Flag Taxa are:</strong></p>",
    "<ul>", paste("<li>", yellowFlag, "</li>", collapse = ""), "</ul>"
  )
  
  HTML(output_text)  # Return the HTML output
})
})
#####THIS IS THE COMP HET PORTION #####
  observeEvent(input$do, {
    req(input$file1)  # Ensure the file input is not NULL
    req(input$comphet)
    req(input$d_type)   # Ensure the data type is selected
        # Read the uploaded FASTA file
    fasta_file <- input$file1$datapath
    data_type <- input$d_type
   
    # Call the function to calculate frequencies and proportions
    frequencies <- calculate_frequency(fasta_file, data_type)

    # Render the outputs for the tables
    output$csrcfv <- renderTable({
      data.frame(Character = names(frequencies$cs_normalized_comparison_sum), 
                 Value = format(frequencies$cs_normalized_comparison_sum, digits=5))
    })
     output$tsrcfv <- renderTable({
      data.frame(Taxa = names(frequencies$ts_normalized_comparison_sum), 
                 Value = format(frequencies$ts_normalized_comparison_sum,digits=5))
    })
    
    output$rcfv <- renderTable({
          data.frame(Character = "RCFV", 
                 Value = format(frequencies$rcfv_normalized_comparison_sum_named,digits=5))
    })
    
    output$cshist <- renderPlot({
       hist(frequencies$cs_normalized_comparison_sum,
            main = "Histogram of ncsRCFV Scores",
            xlab = "csRCFV",
            col = "lightblue",
            border = "black")  # Adjust breaks as necessary
     })

    output$tshist <- renderPlot({
       hist(frequencies$ts_normalized_comparison_sum,
            main = "Histogram of ntsRCFV Scores",
            xlab = "tsRCFV",
            col = "lightblue",
            border = "black")  # Adjust breaks as necessary
     })
    output$ch_bostin <- renderUI({
        tnrcfv_redFlag <- vector()
        tnrcfv_yellowFlag <- vector()
    	ts_mean <- mean(frequencies$ts_normalized_comparison_sum)
    	ts_median <- median(frequencies$ts_normalized_comparison_sum)
    	ts_mad <- mad(frequencies$ts_normalized_comparison_sum)
    	ts_sd <- sd(frequencies$ts_normalized_comparison_sum)
        cnrcfv_redFlag <- vector()
        cnrcfv_yellowFlag <- vector()
    	cs_mean <- mean(frequencies$cs_normalized_comparison_sum)
    	cs_median <- median(frequencies$cs_normalized_comparison_sum)
    	cs_mad <- mad(frequencies$cs_normalized_comparison_sum)
    	cs_sd <- sd(frequencies$cs_normalized_comparison_sum) 
    	for (i in seq_along(frequencies$ts_normalized_comparison_sum)) {
           value <- frequencies$ts_normalized_comparison_sum[i]
           name <- names(frequencies$ts_normalized_comparison_sum)[i]
           if (value > ts_median+(2*ts_mad)) {
                 tnrcfv_yellowFlag <- c(tnrcfv_yellowFlag, name)
                 #print(yellowFlag)
               }
           if (value > ts_median+(2*ts_sd)) {
                 tnrcfv_redFlag <- c(tnrcfv_redFlag, name)
                 #print(redFlag)
               }
            }
        for (i in seq_along(frequencies$cs_normalized_comparison_sum)) {
           value <- frequencies$cs_normalized_comparison_sum[i]
           name <- names(frequencies$cs_normalized_comparison_sum)[i]
           if (value > cs_median+(2*cs_mad)) {
                 cnrcfv_yellowFlag <- c(cnrcfv_yellowFlag, name)
                 #print(yellowFlag)
               }
           if (value > cs_median+(2*cs_sd)) {
                 cnrcfv_redFlag <- c(cnrcfv_redFlag, name)
                 #print(redFlag)
               }
            }
                tagList(
        h4("nRCFV Analysis"),
        p(paste("The total nRCFV of your dataset is:", format(frequencies$rcfv_normalized_comparison_sum_named,digits=5))),
        p("The nRCFV is the average of the taxon-specific nRCFV values and the character-specific nRCFV values."),
        p("This doesn't mean much on its own unless you're comparing between different alignments. If you are, the one with a lower nRCFV has less compositional heterogeneity."),
        p("For further analysis, the histograms show the distribution of tsnRCFV values across the taxa in your dataset and csnRCFV values across the characters in your dataset which, when compared to the average, can help you work out which taxa and characters might end up causing problems in your final phylogenetic analysis. High tsnRCFV and csnRCFV values indicate high levels of compositional heterogeneity."),
        p("We can use the difference between the median and the mean of the tsnRCFV and csnRCFV values to better understand what the distribution of compositional heterogeneity in your dataset looks like."),
        
        h4("Taxon-Specific nRCFV"),
        p("Mean:", format(ts_mean,digits=5)),
        p("Median:", format(ts_median,digits=5)),
        p(ifelse(ts_mean > ts_median,
                 "Your tsnRCFV mean is larger than the median, which indicates that your data contains some sequences that are significantly more compositionally heteogeneous than others. This could be a warning sign.",
                 "Your tsnRCFV mean is smaller than the median, which indicates that your data is distributed towards sequences that are compositionally homogenous. This is a positive sign!")),
        p("From the median, we can identify potentially problematic taxa using two metrics, the standard deviation and the median absolute deviation. In a normal distribution, 95% of data should exist within 2 standard deviations of the mean, but single instances of very large values can inflate the mean by quite a bit. If we instead measure 2 standard deviations from the median, we can more comfortably identify any potential outlier. The median absolute deviation, meanwhile, is a bit like the median's version of a standard deviation. It is the median of how far away each value in the dataset is from the median. It is a stricter measurement, and so all taxa that are 2 standard deviations from the median are marked with a Red Flag so that you can seriously consider excluding them from your analysis, whereas all taxa that are 2 median absolute deviations from the median are marked with a Yellow Flag to warn you that they might be a problem."),
        p("Potential Red Flag taxa are:", paste(tnrcfv_redFlag, collapse = ", ")),
        p("Potential Yellow Flag taxa are:", paste(tnrcfv_yellowFlag, collapse = ", ")),
        
        h4("Character-Specific nRCFV"),
        p("Mean:", format(cs_mean,digits=5)), 
        p("Median:", format(cs_median,digits=5)),
        p(ifelse(cs_mean > cs_median,
                 "Your csnRCFV mean is larger than the median, which indicates that your data contains some sites that are significantly more compositionally heteogeneous than others. This could be a warning sign.",
                 "Your csnRCFV mean is smaller than the median, which indicates that your data is distributed towards sites that are compositionally homogenous. This is a positive sign!")),
         p("Like with the tsnRCFV, we can assess the standard deviation and median absolute deviation from the median to establish whether there are any problematic characters. Unfortunately, how to deal with difficult characters is a bit more complex!"),        
         p("Potential Red Flag characters are:", paste(cnrcfv_redFlag, collapse = ", ")),
         p("Potential Yellow Flag characters are:", paste(cnrcfv_yellowFlag, collapse = ", ")),
         if (length(cnrcfv_yellowFlag) > 4 & data_type=="AA"){
	          p("As 5 or more characters show high compositional heterogeneity, it might be a good idea to consider an approach such as Dayhoff 6-state recoding, and comparing the trees produced by that process against ones produced by your data without any modification. Be aware, while simplifying the amino acid alphabet down from 20 to 6 characters reduces compositional heterogeneity, it might also mask useful phylogenetic information! Check Hernandez et al 2019")
         },
        if (length(cnrcfv_redFlag) > 1 & data_type=="DNA"){
	          p("As 2 or more characters show high compositional heterogeneity, it might be a good idea to consider an approach such as 2-state recoding to purines a pyrimidines, and comparing the trees produced by that process against ones produced by your data without any modification. It might be the result of a significant AT or CG bias in your data, which is sometimes biological, and sometimes artefactual. Be aware, while simplifying the nucleotide alphabet down from 4 to 2 characters reduces compositional heterogeneity, it might also mask useful phylogenetic information! Check Hernandez et al 2019")
         },
         if (length(cnrcfv_redFlag) > 0){
            	p("Here, only a few characters show high levels of compositional heterogeneity. It might be worth considering a site-sensitive model, such as CAT+GTR for your analysis. You might not need to employ a recoding strategy, as it is possible that you will lose information by masking the variation within your recoded groups.")
         }
         else{
             	p("There weren't any particularly problematic characters identified. It doesn't seem like site to site compositional heterogeneity is going to be a problem for your analysis!")
         }

)})
})

######THIS IS THE SITE SATURATION PORTION########
  observeEvent(input$do, {
    req(input$file1)
    req(input$sitesat)
    req(input$d_type)
    data_type <- input$d_type
    ######THIS IS THE DE SCORE PORTION#####
    if(data_type=="AA"){
    # Read the uploaded FASTA file
    fasta_file <- input$file1$datapath
    de_score <- py$calculate_deScore(fasta_file)
#    print(de_score$TotalFreqResults)
	    # Render the tables in the UI
    total_freq_results <- as.data.frame(de_score$TotalFreqResults, stringsAsFactors = FALSE)
    
    # Convert TaxResults (Python dict) to a data frame
    tax_results_list <- lapply(de_score$TaxResults, function(x) {
      # Convert each tax result to a data frame
      as.data.frame(x, stringsAsFactors = FALSE)
    })
    
    # Combine tax result data frames into one
    tax_results_df <- do.call(rbind, tax_results_list)
    tax_scores <- sapply(de_score$TaxResults, function(x) x$DE_Score)
        
    # Render the tables in the UI
    output$freq_table <- renderTable({
      total_freq_results
    })

    output$sathist <- renderPlot({
       hist(tax_scores,
            main = "Histogram of tDE-Scores",
            xlab = "taxon-specific DE Score",
            col = "lightblue",
            border = "black")  # Adjust breaks as necessary
     })

    output$taxa_table <- renderTable({
      tax_results_df
    })

  output$sitesat_bostin <- renderUI({
        tde_redFlag <- vector()
        tde_yellowFlag <- vector()
    	tde_mean <- mean(tax_results_df$DE_Score)
    	tde_median <- median(tax_results_df$DE_Score)
    	for (i in seq_along(tax_results_df$DE_Score)) {
           value <- tax_results_df$DE_Score[i]
           name <-  tax_results_df$Sequence[i]
           if (value < 0.355) {
                 tde_yellowFlag <- c(tde_yellowFlag, name)
                 #print(yellowFlag)
               }
           if (value < 0) {
                 tde_redFlag <- c(tde_redFlag, name)
                 #print(redFlag)
               } 
            }
                tagList(
        h4("DE-Score Analysis"),
        p(paste("The DE-Score of your dataset is:", format(total_freq_results$DE_Score,digits=5))),
        p("The DE-Score is the average of the taxon-specific DE-Scores."),
        p("When evaluating the total DE-Score, consider two things. First, a DE-Score below 0 indicates that the data is completely saturated, and likely uninformative - for best results, seek a DE-Score of 0.354 or higher. This means that the alignment is 2 minimum-information steps away from total saturaton (2x0.177 away from DE-Score 0). When comparing between alignments, the one with a higher DE-Score has less entropic site saturation."),
        p("For further analysis, the histograms show the distribution of taxon-specific DE-Scores across the taxa in your dataset, which can help you work out which taxa might end up causing problems in your final phylogenetic analysis."),
        p("We can use the median and the mean of the t-DEScore to better understand what the distribution of site saturation in your dataset looks like."),
        
        h4("Taxon-Specific DE-Scores"),
        p("Mean:", format(tde_mean,digits=5)),
        p("Median:", format(tde_median,digits=5)),
        p(ifelse(tde_mean < tde_median,
                 "Your tDE-Score mean is smaller than the median, which indicates that your data may contain some sequences that are significantly more saturated than others. This could be a warning sign.",
                 "Your tDE-Score mean is larger than the median, which indicates that your data is distributed towards sequences that are unsaturated. This is a positive sign!")),
        p("In addition to this, we can identify problematic taxa by observing their distance from the DE-Score's critical value of 0. tDE-Scores below 0 are marked as red flag taxa, and might be best to remove, whereas those below 0.354 - two uninformative steps away from the critical value - are treated as yellow flags, and should be considered with caution."),
        p("Potential Red Flag taxa are:", paste(tde_redFlag, collapse = ", ")),
        p("Potential Yellow Flag taxa are:", paste(tde_yellowFlag, collapse = ", ")),
)}
)
}
	else{
	        # File input handling
        fasta_file <- input$file1$datapath
        
        # Call the Python function using reticulate
        c_score <- py$calculate_CScore(fasta_file)
        
        total_freq_results <- as.data.frame(c_score$TotalFreqResults, stringsAsFactors = FALSE)
    
        # Convert TaxResults (Python dict) to a data frame
        tax_results_list <- lapply(c_score$TaxResults, function(x) {
        # Convert each tax result to a data frame
            as.data.frame(x, stringsAsFactors = FALSE)
        })
    
        # Combine tax result data frames into one
        tax_results_df <- do.call(rbind, tax_results_list)
        tax_scores <- sapply(c_score$TaxResults, function(x) x$CScore)
        # Prepare output for frequency table
        output$freq_table <- renderTable({
            total_freq_results
        })

       output$sathist <- renderPlot({
       hist(tax_scores,
            main = "Histogram of C Values",
            xlab = "C Value",
            col = "lightblue",
            border = "black")  # Adjust breaks as necessary
     })

        # Prepare output for taxa table
    output$taxa_table <- renderTable({
       tax_results_df
    })
    output$sitesat_bostin <- renderUI({
        tC_redFlag <- vector()
        tC_yellowFlag <- vector()
    	tC_mean <- mean(tax_results_df$CScore)
    	tC_median <- median(tax_results_df$CScore)
    	for (i in seq_along(tax_results_df$CScore)) {
           value <- tax_results_df$CScore[i]
           name <-  tax_results_df$Sequence[i]
           if (value < 20) {
                 tC_yellowFlag <- c(tC_yellowFlag, name)
                 #print(yellowFlag)
               }
           if (value < 10) {
                 tC_redFlag <- c(tC_redFlag, name)
                 #print(redFlag)
               } 
            }
                tagList(
        h4("C Factor"),
	    p("As you are using amino acid data, we selected the C Factor - the Convergence Factor. This is based on the ratio of the standard deviation of the Transition/Transversion ratio of the dataset and the standard deviation of the uncorrected p-distance. It also has the advantage that it can be used to detect taxa that are particularly saturated. We call this the taxon C Factor."),
	    p(paste("The C Factor of this entire dataset is", format(total_freq_results$CScore, digits=5), ". Meanwhile, the standard deviation of the transition transversion ratios is", format(total_freq_results$TiTv.StD, digits=5), "and the standard deviation of the uncorrected p-distances is", format(total_freq_results$PValue, digits=5), ".")),
		if (total_freq_results$CScore > 20) {
           p("As the C Factor is greater than 20, this means that saturation is unlikely to be a problem with this dataset, as the Observed and Expected Ti/Tv ratio is likely greater than 1 (Struck et al 2008).")
           }
        else if(total_freq_results$CScore > 10){
	       p("As the C Factor is between 10 and 20, it means that this dataset is not entirely saturated, but that saturation might affect topological reconstruction. As it is still quite low, you might want to observe the distribution of the taxa C Factors, to see if any taxa are contributing particularly to the overall score.")
           }
        else{
	       p("As the C Factor is below 10, this means that saturation might be a serious problem in this dataset. Check the taxa C Factors to see which taxa are contributing to this.")
           },
        p("As with the DE-Score, we can use the taxon-specific C Factors to directly assess individual taxa for entropic site saturation in your dataset. Sequences with a taxon C Factor of 20 or lower are marked as yellow flag taxa, and should be treated with caution. Sequences with a taxon C Factor of 10 or lower are instead marked as red flag taxa, and may cause problems with phylogenetic reconstruction."),
        p("Potential Red Flag taxa are:", paste(tC_redFlag, collapse = ", ")),
        p("Potential Yellow Flag taxa are:", paste(tC_yellowFlag, collapse = ", ")),

)})}
})
}