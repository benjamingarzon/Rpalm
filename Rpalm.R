# An R wrapper for FSL's palm (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/)
# benjamin.garzon@gmail.com
# February 2018

PALM_DIR='~/Software/palm-alpha109/'

Rpalm = function(mydata, dependent, myformula, mycontrasts, options = '', output_dir = NULL, kinship_file  = NULL, contrast_type = 'both'){
  
  mywd = getwd()
  # parse the formula
  all = strsplit(myformula, '~')[[1]]
  independent = strTrim(unlist(strsplit(all[2], '\\+')))
  if (!is.null(kinship_file)) mydata = subset(mydata, ZygositySR !=' ')
  mydata = mydata[complete.cases(mydata[c("Subject", dependent, independent)]), ]
  mydata = as.data.frame(lapply(mydata, as.numeric))
  Y = mydata[, dependent]
  X = as.data.frame(scale(mydata[, independent]))

  X$Intercept = 1
  if (is.null(output_dir)) {
    output_dir = tempdir()
    cleanafter = T
    
  } else {
    if (file.exists(output_dir)) unlink(output_dir, recursive = T)
    cleanfter = F
  } 
  
  dir.create(output_dir)
  print(paste("Creating new directory: ", output_dir))
  setwd(output_dir)
  
  numwaves = ncol(X)
  numpoints = nrow(X)
  numcontrasts = length(mycontrasts)

  # data 
  write.table(Y, file = 'data.csv', quote = F, row.names = F, col.names = F, sep = ',')
  
  X = X[, c(numwaves, seq(numwaves - 1))]
  
  # design matrix
  cat(
    '/NumWaves',  '\t', numwaves,  '\n', 
    '/NumPoints', '\t', numpoints, '\n', 
    '\n',  
    '/Matrix\n', 
    file='design.mat', sep='')
  
  write.table(X, file = 'design.mat', quote = F, row.names = F, col.names = F, append = T, sep = '\t')
  
  # contrasts
  numcontrasts = length(mycontrasts)
  contrast.matrix = matrix(0, numcontrasts, numwaves)
  colnames(contrast.matrix) = colnames(X)
  contrast.matrix[ , mycontrasts] = eye(numcontrasts)

    
  if (contrast_type == 'both') {
    contrast.matrix = rbind(contrast.matrix, -contrast.matrix)
    numcontrasts = nrow(contrast.matrix)
    
    rownames(contrast.matrix) = c(paste0(mycontrasts, ">0"), paste0(mycontrasts, "<0"))
    # reorder the matrix
    contrast.matrix = rbind(
      contrast.matrix[seq(1, numcontrasts, 2), ],
      contrast.matrix[seq(2, numcontrasts, 2), ]
    )
    
  }

  if (contrast_type == '<') {
    contrast.matrix = -contrast.matrix
    rownames(contrast.matrix) = paste0(mycontrasts, "<0")
  }
  
  if (contrast_type == '>') {
    rownames(contrast.matrix) = paste0(mycontrasts, ">0")
  }
  
  numcontrasts = nrow(contrast.matrix)
  
  colnames(contrast.matrix) = colnames(X)

  
  print("Contrast matrix: ")
  print(contrast.matrix)
  
  cat(
    paste(paste(paste0('/ContrastName', seq(numcontrasts)), 
                rownames(contrast.matrix), sep = '\t'), collapse = '\n'),   
    '\n',  
    '/NumWaves',  '\t', numwaves,  '\n', 
    '/NumContrasts', '\t', numcontrasts, '\n', 
    '\n',  
    '/Matrix\n', 
    file='design.con', sep='')
  
  write.table(contrast.matrix, file = 'design.con', quote = F, row.names = F, col.names = F, append = T, sep = '\t')
  
  # generate blocks file
  if (!is.null(kinship_file)){
    system(
      paste0('matlab -nodisplay -nosplash -nodesktop -r ',
             '"addpath ', PALM_DIR, ';',
             'hcp2blocks(\'', kinship_file, 
             '\',\'EB.csv\', true, [', paste(mydata$Subject, collapse = ' '), ']\'); exit;"')
    ) 
    options = paste(options, '-eb EB.csv')
  }

  # run palm  
  PALM = file.path(PALM_DIR, 'palm')
  system('rm results_dat_tstat*')
  command = paste(PALM, '-i data.csv -d design.mat -t design.con -o results', options, ' > log')
  print(command)
  system(command)
  
  # read out results
  if (numcontrasts > 1) {
    system('cat results_dat_tstat_c*.csv > results_dat_tstat.csv')
    system('cat results_dat_tstat_uncp_c*.csv > results_dat_tstat_uncp.csv')
    system('cat results_dat_tstat_fwep_c*.csv > results_dat_tstat_fwep.csv')
  }
  tstat = read.table('results_dat_tstat.csv', sep =',')
  uncp = read.table('results_dat_tstat_uncp.csv', sep =',')
  fwep = read.table('results_dat_tstat_fwep.csv', sep =',')
  colnames(tstat) = colnames(uncp) = colnames(fwep) = dependent
  rownames(tstat) = rownames(uncp) = rownames(fwep) = rownames(contrast.matrix)
  
  print('FWE-corrected p-values')
  print(fwep)
  # clean up
  if (F) {
    setwd(mywd)
    unlink(output_dir, recursive = T)
  }
  results = list(tstat = tstat,
                 uncp = uncp,
                 fwep = fwep
                 )
  return(results)
}

