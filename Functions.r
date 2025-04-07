calculate_LBi <- function(input_fasta){
      DistMatrix <- dist.ml(input_fasta)
      treeNJ <- NJ(DistMatrix)
      dist.mat <- cophenetic.phylo(treeNJ)
      dist.mat[dist.mat == 0] <- NA
      matrixmean <- mean(dist.mat, na.rm = TRUE)
      result_matrix <- matrix(NA, nrow = nrow(dist.mat), ncol = 2)
      
      # Calculating average pairwise distance and LB-Score
      for (col in 1:ncol(dist.mat)) {
        t_avgdist <- mean(dist.mat[col,], na.rm = TRUE)
        t_deviation <- ((t_avgdist / matrixmean) - 1) * 100
        result_matrix[col, ] <- c(t_avgdist, t_deviation)
      }
      
      rownames(result_matrix) <- rownames(dist.mat)
      colnames(result_matrix) <- c("averagePDist", "LB-Score")
      return(list(result_matrix = result_matrix, dist.mat = dist.mat, treeNJ = treeNJ))
}

calculate_frequency <- function(fasta_file, data_type) {
  # Read the FASTA file using Biostrings package
  if (data_type == "DNA") {
    sequences <- readDNAStringSet(fasta_file)
  } else if (data_type == "AA") {
    sequences <- readAAStringSet(fasta_file)
  } else {
    stop("Invalid data type. Please choose 'nucleotides' or 'amino acids'.")
  }
  
  # Initialize an empty list to store the frequency counts
  freq_table <- list()
  taxa_names <- names(sequences)
  valid_nucleotides <- c('A', 'T', 'C', 'G')  # Include gap character for aligned sequences
  valid_amino_acids <- c('A','V','L','I','P','M','F','W','G','S','T','C','N','Q','Y','D','E','K','R','H')  # 20 amino acids + gap
  proportion_table <- list()
  total_length <- 0
  total_sequences <- length(sequences)  # Total number of sequences
  alignment_length <- nchar(as.character(sequences[1]))  # Length of the alignment (length of a single sequence)

  # Define the valid characters for nucleotides and amino acids
 
  # Clean sequences (remove invalid characters)
  if (data_type == "DNA") {
     overall_freqs <- rep(0, length(valid_nucleotides))  # For nucleotide sequences
     sequences <- gsub(paste0("[^", paste(valid_nucleotides, collapse = ""), "]"), "", as.character(sequences))  # For nucleotide sequences
     } else if (data_type == "AA") {
     overall_freqs <- rep(0, length(valid_amino_acids))
     sequences <- gsub(paste0("[^", paste(valid_amino_acids, collapse = ""), "]"), "", as.character(sequences))  # For amino acids
  }
  
    
  # Iterate over each sequence (now represented as character strings)
  for (i in seq_along(sequences)) {
    seq_str <- as.character(sequences[i])
    taxa_name <- taxa_names[i]
    
    # Get the length of the sequence (including gaps)
    seq_length <- nchar(seq_str)
    total_length <- total_length + seq_length
    
    # Count frequencies of characters (nucleotides or amino acids)
    if (data_type == "DNA") {
      freqs <- table(factor(strsplit(seq_str, NULL)[[1]], levels = valid_nucleotides))
      overall_freqs <- overall_freqs + freqs
    } else if (data_type == "AA") {
      freqs <- table(factor(strsplit(seq_str, NULL)[[1]], levels = valid_amino_acids))
      overall_freqs <- overall_freqs + freqs
    }
    
    # Convert frequencies to proportions (divide by sequence length)
    proportions <- freqs / seq_length
    
    # Store the proportions in the table with the taxa name
    proportion_table[[taxa_name]] <- proportions
  }
  overall_proportions <- overall_freqs / total_length
    # Compare each taxa's proportions with the dataset's overall proportions
  comparison_table <- list()
  
  for (taxa_name in names(proportion_table)) {
    taxa_proportions <- proportion_table[[taxa_name]]
    
    # Calculate the absolute difference between taxa proportions and overall proportions
    comparison <- abs(taxa_proportions - overall_proportions)
    
    # Divide each value by the total number of sequences
    comparison_normalized <- comparison / total_sequences
    
    comparison_table[[taxa_name]] <- comparison_normalized
  }
  # Now normalize the total_comparison_sum by the requested formula
  # Compute the normalization factor
  num_taxa <- total_sequences  # Number of taxa (sequences)
  num_characters <- if (data_type == "DNA") length(valid_nucleotides) else length(valid_amino_acids)  # Number of characters (nucleotides or amino acids)
  
  # Normalization factor for RCFV
  rcfv_normalization_factor <- (alignment_length ^ -0.5) * (num_taxa ^ 0.01) * num_characters * 100
  cs_normalization_factor <- (alignment_length ^ -0.5) * 100
  ts_normalization_factor <- (alignment_length ^ -0.5) * (num_taxa ^ -1) * num_characters * 100

  # 1. Sum of all comparison_normalized for each character (nucleotide/amino acid)
  comparison_sums_by_character <- sapply(names(comparison_table[[1]]), function(char) {
    sum(sapply(comparison_table, function(comp) comp[char]))
  })
  cs_normalized_comparison_sum <- comparison_sums_by_character / cs_normalization_factor
  
  # 2. Sum of all comparison_normalized for each sequence (taxon)
  comparison_sums_by_sequence <- sapply(names(comparison_table), function(taxa_name) {
    sum(comparison_table[[taxa_name]])
  })
  ts_normalized_comparison_sum <- comparison_sums_by_sequence / ts_normalization_factor

  # 3. Sum of all comparison_normalized for the entire dataset
  total_comparison_sum <- sum(unlist(comparison_sums_by_sequence))
  
  # Normalize the total_comparison_sum
  rcfv_normalized_comparison_sum_named <- total_comparison_sum / rcfv_normalization_factor
#  result_name <- "nRCFV"  # You can specify whatever name you want here
#  rcfv_normalized_comparison_sum_named <- paste(result_name, ":", rcfv_normalized_comparison_sum)
  
  # Return both the proportion table, comparison table, and sums
  return(list(cs_normalized_comparison_sum = cs_normalized_comparison_sum,
              ts_normalized_comparison_sum = ts_normalized_comparison_sum,
              rcfv_normalized_comparison_sum_named = rcfv_normalized_comparison_sum_named))
}

py_run_string("from Bio import SeqIO
import numpy as np

def calculate_deScore(fasta_file):
    # Read the sequences using Biopython's SeqIO
    alignFile = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    
    # Amino acid groups
    Small = set(['A', 'G', 'P', 'S', 'T'])
    AcidAmide = set(['D', 'E', 'N', 'Q'])
    Basic = set(['H', 'K', 'R'])
    Hydrophobic = set(['I', 'L', 'V', 'M'])
    Aromatic = set(['F', 'W', 'Y'])
    Sulfur = set(['C'])
    
    # Initialize lists for storing results
    PhySeqs = {}
    AllTvs = []
    AllTis = []
    AllTiFreqs = []
    TaxResults = {}
    
    # Convert sequences into a list where each element is a sequence
    for seqID, seqRecord in alignFile.items():
        PhySeqs[seqID] = list(str(seqRecord.seq))
    
    # Get the size of the dataset
    size = len(PhySeqs)
    
    # Loop over each sequence (k) and compare it with all other sequences (compk)
    for k, vals in PhySeqs.items():
        TaxTvs = []
        TaxTis = []
        TaxTiFreqs = []
        
        for compk, compvals in PhySeqs.items():
            if k == compk: continue  # Skip comparing the sequence to itself
            
            pairwiseTI = 0
            pairwiseTV = 0
            SameCounter = 0
            GapCounter = 0
            
            # Compare the sequences site by site
            for i, (query, compquery) in enumerate(zip(vals, compvals)):
                if query == '-' or compquery == '-':
                    GapCounter += 1
                elif query == compquery:
                    SameCounter += 1
                elif query in Small and compquery in Small:
                    pairwiseTI += 1
                elif query in AcidAmide and compquery in AcidAmide:
                    pairwiseTI += 1
                elif query in Basic and compquery in Basic:
                    pairwiseTI += 1
                elif query in Hydrophobic and compquery in Hydrophobic:
                    pairwiseTI += 1
                elif query in Aromatic and compquery in Aromatic:
                    pairwiseTI += 1
                elif query in Sulfur and compquery in Sulfur:
                    pairwiseTI += 1
                else:
                    pairwiseTV += 1
            
            # Calculate TiFrequencies and store values
            if pairwiseTI + pairwiseTV > 0:
                TiFreq = pairwiseTI / (pairwiseTI + pairwiseTV)
                TaxTis.append(pairwiseTI)
                TaxTvs.append(pairwiseTV)
                TaxTiFreqs.append(TiFreq)
                AllTis.append(pairwiseTI)
                AllTvs.append(pairwiseTV)
                AllTiFreqs.append(TiFreq)
        
        # Calculate average and standard deviation for taxa
        TaxAverageTransI = np.mean(TaxTis) if TaxTis else 0
        TaxAverageTransV = np.mean(TaxTvs) if TaxTvs else 0
        TaxAverageTiFreq = np.mean(TaxTiFreqs) if TaxTiFreqs else 0
        TaxStdTransI = np.std(TaxTis) if TaxTis else 0
        TaxStdTransV = np.std(TaxTvs) if TaxTvs else 0
        TaxStdTiFreq = np.std(TaxTiFreqs) if TaxTiFreqs else 0
        TaxDist = TaxAverageTiFreq - 0.177
        TaxDE = TaxDist / (0.255 * (size ** -0.15))
        
        # Store results in TaxResults
        TaxResults[k] = {
            'Sequence': k,
            'AverageTiFreq': TaxAverageTiFreq,
            'StdTiFreq': TaxStdTiFreq,
            'DE_Score': TaxDE
        }
    
    # Calculate the overall averages and standard deviations
    AverageTransI = np.mean(AllTis)
    AverageTransV = np.mean(AllTvs)
    AverageTiFreq = np.mean(AllTiFreqs)
    StdTransI = np.std(AllTis)
    StdTransV = np.std(AllTvs)
    StdTiFreq = np.std(AllTiFreqs)
    FreqDist = AverageTiFreq - 0.177
    totalDE = FreqDist / (0.255 * size ** -0.15)
    
    # Store total frequency results
    TotalFreqResults = {
        'FileName': 'DE-Score',
        'AverageTiFreq': AverageTiFreq,
        'StdTiFreq': StdTiFreq,
        'DE_Score': totalDE
    }
    
    return {'TotalFreqResults': TotalFreqResults, 'TaxResults': TaxResults}"
)

py_run_string("
import numpy as np
from Bio import SeqIO

def calculate_CScore(fasta_file):
    # Read the sequences using Biopython's SeqIO
    alignFile = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
        
    # Define purines and pyrimidines
    Purines = ['A', 'G']
    Pyrimidines = ['T', 'C']
    
    # Initialize lists for storing results
    PhySeqs = {}
    AllTiTvs = []
    AllPs = []
    AllCs = []
    taxa_data = {}
    TaxResults = {}
        
    # Convert sequences into a list where each element is a sequence
    for seqID, seqRecord in alignFile.items():
        PhySeqs[seqID] = list(str(seqRecord.seq))
        
    # Get the size of the dataset
    size = len(PhySeqs)
        
    # Loop over each sequence (k) and compare it with all other sequences (compk)
    for k, vals in PhySeqs.items():
        TaxTvs = []
        TaxTis = []
        TaxTiTvs = []
        TaxPs = []
        TaxCs = []
            
        for compk, compvals in PhySeqs.items():
            if k == compk: continue  # Skip comparing the sequence to itself
                
            pairwiseTI = 0
            pairwiseTV = 0
            SameCounter = 0
            GapCounter = 0
                
            # Compare the sequences site by site
            for i, (query, compquery) in enumerate(zip(vals, compvals)):
                if query == '-' or compquery == '-':
                    GapCounter += 1
                elif query == compquery:
                    SameCounter += 1
                elif query in Purines and compquery in Purines:
                    pairwiseTI += 1
                elif query in Pyrimidines and compquery in Pyrimidines:
                    pairwiseTI += 1
                elif query in Purines and compquery in Pyrimidines:
                    pairwiseTV += 1
                elif query in Pyrimidines and compquery in Purines:
                    pairwiseTV += 1
                else:
                    AmbigCounter += 1
            # Calculate TiFrequencies and store values
            if pairwiseTI == 0 or pairwiseTV == 0:
            	TaxTiTv = 0
            	Taxp_value = 0
            	Taxc_value = 0
            else:
                TaxTiTv = pairwiseTI / pairwiseTV
                Taxp_value = (pairwiseTI + pairwiseTV)/(pairwiseTI + pairwiseTV + SameCounter)
                TaxTiTvs.append(TaxTiTv)
                TaxPs.append(Taxp_value)
                AllTiTvs.append(TaxTiTv)
                AllPs.append(Taxp_value)

        
        # Calculate average and standard deviation for taxa
        TaxStdTiTv = np.std(TaxTiTvs) if TaxTiTvs else 0
        TaxStdPs = np.std(TaxPs) if TaxPs else 0
        TaxStdCs = TaxStdTiTv/TaxStdPs
        
        # Store results in TaxResults
        TaxResults[k] = {
            'Sequence': k,
            'TiTv_Std': TaxStdTiTv,
            'PValue': TaxStdPs,
            'CScore': TaxStdCs
        }
    
    # Calculate the overall averages and standard deviations
    StdTiTv = np.std(AllTiTvs)
    StdP = np.std(AllPs)
    StdC = StdTiTv/StdP
    
    # Store total frequency results
    TotalFreqResults = {
            'FileName': 'C Score',
            'TiTv StD': StdTiTv,
            'PValue': StdP,
            'CScore': StdC
        }
    
    return {'TotalFreqResults': TotalFreqResults, 'TaxResults': TaxResults}"
)