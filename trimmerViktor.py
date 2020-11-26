## Import module with functions required for the program
import trimmer_functions as tf
import gzip
import sys
import os











## Get filenames and arguments from command line
args = tf.run_arg_parser()

try:
    LEADING = int(args.LEADING)               # 3' bases to be removed                                      ## ERROR if user input is unexpected
    TRAILING = int(args.TRAILING)             # 5' bases to be removed
    BASE_QUALITY = int(args.BASEQUALITY)      # Quality threshold for single bases (3 by default)
    AVG_QUALITY = int(args.AVGQUALITY)        # Average quality threshold for read/window (15 by default)
    MIN_LEN = int(args.MINLEN)                # Minimum length for trimmed read (36 by default)
    N_MAX = int(args.MAXN)                    # Maximum number of unknown bases in read (3 by default)      
    WIN_SIZE = int(args.WINDOWSIZE)           # Window size for sliding window approach (1 by default)
    USER_PHRED = args.PHRED                   # User given phred type                                     
except ValueError as err:
    print('Invalid input. Reason: ' + str(err))
    sys.exit(1)

in_fwFile = args.FILE1
in_revFile = args.FILE2

if USER_PHRED not in ['', '33', '64']:                                                                  ##### DEAL WITH 33 vs phred33 vs phred 33 vs PHRED33    
    print('Invalid input for phred type: {} \nAccepted input: \'33\', \'64\''.format(USER_PHRED))
    sys.exit(1)

# Filename base
base_fw = in_fwFile.split('.')[0]
base_rev = in_revFile.split('.')[0]


###################
# Singe end reads #
################### 
if in_revFile == '':
    print('Single end!')
    
    ## -------- Trimmer ---------- ##
    try:
        # Determine Phred encoding type
        phred = tf.phred_control(in_fwFile, USER_PHRED) 
    
        # Open files, checking if file is compressed or not 
        if in_fwFile.endswith('.gz'):                                                   ## ERROR with file opening (read and write - if input files don't exist or if output files already exist (not ot overwrite them))
            file_fw = gzip.open(in_fwFile, 'rt')
        else:    
            file_fw = open(in_fwFile, 'r')
    except IOError as err:
        print('File could not be opened. Reason: ' + str(err))
        sys.exit(1)
    
    # The output file, controlling if file already exists in woriking directory
    out_fw = tf.controling_output_file(in_fwFile)                                       ### Added this function to prevent overwriting files
    
    # Read first line and initialize variables
    line_fw = file_fw.readline()
    
    line_count = 0
    fastq_fw = []

    # Variables for stats
    dropped_reads = 0
    trimmed_reads = 0
    if (LEADING != 0) or (TRAILING != 0):
        trimmed_reads = 'all'    
    read_count = 0
    read_len_sum = 0
    read_qual_sum = 0
    trimmed_read_len_sum = 0
    trimmed_read_qual_sum = 0

    # Iterate through lines
    while line_fw != '':
        line_count += 1
        
        # For each read, make a 4 element list with each line as an element
        fastq_fw.append(line_fw.strip())

        if line_count == 4:
            # Keep track of number of reads in files                        #### ERROR if number of reads are not the same (corrupt files) -- in the end of the loop, with the .readline()
            read_count += 1

            # Read sequence and quality (encoded and decoded)
            read_fw = fastq_fw[1]
            qual_str_fw = fastq_fw[3]
            qual_score_fw = tf.quality_score(qual_str_fw, phred)

            ## Stats for original reads ##
            # Average length of all kept reads
            read_len_sum += len(read_fw)
            # Average quality of all kept reads
            avg_qual = sum(qual_score_fw)/len(qual_score_fw)
            read_qual_sum += avg_qual          

            ## STEP 0: Drop read if quality can't be determined            
            if qual_score_fw == 'unknown':
                # Keep track of dropped read
                dropped_reads += 1
                # Reset
                line_count, fastq_fw = 0, []
                # Read new line
                line_fw = file_fw.readline()
                # Continue to next read without printing                
                continue
            
            ## STEP 1: Remove leading and trailing bases, given user input
            read_fw, qual_str_fw, qual_score_fw = tf.removal_of_bases(read_fw, qual_str_fw, qual_score_fw, LEADING, TRAILING)
            
            ## STEP 2: Remove leading and trailing bases, based on quality
            read_fw, qual_str_fw, qual_score_fw, trimmed = tf.sliding_window_pop(read_fw, qual_str_fw, qual_score_fw, WIN_SIZE, AVG_QUALITY, BASE_QUALITY)
            if (trimmed == True) and (trimmed_reads != 'all'):
                trimmed_reads += 1

            ## STEP 3: Drop reads that become too short after trimming
            if (len(read_fw) < MIN_LEN):                                    
                # Keep track of dropped read
                dropped_reads += 1
                # Reset
                line_count, fastq_fw = 0, []
                # Read new line
                line_fw = file_fw.readline()
                # Continue to next read without printing
                continue

            ## STEP 4: Drop reads with low average quality
            avg_qual_fw = sum(qual_score_fw)/len(qual_score_fw)
            if (avg_qual_fw < AVG_QUALITY):
                # Keep track of dropped read
                dropped_reads += 1
                # Reset
                line_count, fastq_fw = 0, []
                # Read new line
                line_fw = file_fw.readline()
                # Continue to next read without printing
                continue

            # STEP 5: Drop reads with too many N bases
            if (read_fw.count('N') > N_MAX):
                # Keep track of dropped read
                dropped_reads += 1                
                # Reset
                line_count, fastqForward = 0, []
                # Read new lines
                line_fw = file_fw.readline()    
                # Continue to next read without printing
                continue 

            # STEP 5: Print trimmed reads onto outfile
            tf.print_read(fastq_fw[0], read_fw, qual_str_fw, out_fw)

            ## Stats for trimmed reads ##
            # Average length of all kept reads
            trimmed_read_len_sum += len(read_fw)
            # Average quality of all kept reads
            trimmed_read_qual_sum += avg_qual_fw           

            # Reset
            line_count, fastq_fw = 0, []

        # Read next line
        line_fw = file_fw.readline()

    # Close files
    file_fw.close()
    out_fw.close()
    
    ## --------------------------- ##

    ## ------- Log file ---------- ##
    log = open(base_fw + '.log', 'w')

    print('This is the log file for the trimming of', in_fwFile, file=log)

    print('\n===============\nSETTINGS\n===============', file=log)
    print('File 1:', in_fwFile, file=log)                 # File 1                                      
    print('Base quality:', BASE_QUALITY, file=log)        # Single base quality
    print('Average quality:', AVG_QUALITY, file=log)      # Average quality
    print('Lead trim:', LEADING, file=log)                # Lead trim
    print('Trail trim:', TRAILING, file=log)              # Trail trim
    print('Window size:', WIN_SIZE, file=log)             # Window size 
    print('Maximum unknown bases:', N_MAX, file=log)      # Maximum number of Ns
    print('Min lenght:', MIN_LEN, file=log)               # Minimum lenght after trim
    if (USER_PHRED != '') and (phred != USER_PHRED):      # Phred encoding type
        print("Phred encoding was set to {}, but {} was used.".format(USER_PHRED, phred), file=log)
    else:
        print('Phred: ' + phred, file=log)            

    print('\n===============\nSTATS\n===============', file=log)
    print('*** Before trimming ***', file=log)
    print('Total number of reads:', read_count, file=log)
    print('Average length of reads:', round(read_len_sum/read_count,2), file=log)
    print('Average quality of reads:', round(read_qual_sum/read_count,2), file=log)
    print('\n*** After trimming ***', file=log)
    print('Read pairs dropped due to low quality/short length:', dropped_reads, file=log)
    print('Read pairs kept:', read_count-dropped_reads, file=log)
    print('Reads trimmed:', trimmed_reads, file=log)
    print('Average length of kept reads:', round(trimmed_read_len_sum/(read_count-dropped_reads),2), file=log)
    print('Average quality of kept reads:', round(trimmed_read_qual_sum/(read_count-dropped_reads),2), file=log)

    log.close()    

    ## --------------------------- ##

    ## -------- STDOUT --------- ##
    # If trimming was successful
    print('Congratulations! Your trimming was successful. \nYou can find your results in the', base_fw + '_trimmed.fa file and some additional info in the the', base_fw + '.log file. \nPleasure working with you!')
    ## ------------------------- ##



####################
# Paired end reads #
####################
else:
    print('Paired end!')
    
    ## -------- Trimmer ---------- ##
    try:
        # Check whether both files have the same encoding
        if tf.phred_control(in_fwFile, USER_PHRED) != tf.phred_control(in_revFile, USER_PHRED):
            print("ERROR, the two given files do not have the same encoding type. ")
            sys.exit(1)

        # Determine Phred encoding type
        phred = tf.phred_control(in_fwFile, USER_PHRED)

        # Open files, checking if file is compressed or not 
        if in_fwFile.endswith('.gz'):
            file_fw = gzip.open(in_fwFile, 'rt')
        else:    
            file_fw = open(in_fwFile, 'r')

        if in_revFile.endswith('.gz'):
            file_rev = gzip.open(in_revFile, 'rt')
        else:    
            file_rev = open(in_revFile, 'r')
    except IOError as err:
        print('File could not be opened. Reason: ' + str(err))
        sys.exit(1)

    # Output files, controlling if file already exists in working directory                     ########## Added this function
    out_fw = tf.controling_output_file(in_fwFile)
    out_rev = tf.controling_output_file(in_revFile)

    # Read first line of each file and initialize variables
    line_fw = file_fw.readline()
    line_rev = file_rev.readline() 

    line_count, fastq_fw, fastq_rev = 0, [], []

    # Variables for stats
    dropped_reads = 0
    trimmed_reads = 0
    if (LEADING != 0) or (TRAILING != 0):
        trimmed_reads = 'all'    
    read_count = 0
    read_len_sum = 0
    read_qual_sum = 0
    trimmed_read_len_sum = 0
    trimmed_read_qual_sum = 0

    # Iterate through lines
    while (line_fw or line_rev) != '':
        line_count += 1
        
        # For each file, and each read, make a 4 element list with each line as an element
        fastq_fw.append(line_fw.strip())
        fastq_rev.append(line_rev.strip())

        if line_count == 4:
        
            # Keep track of number of reads in files                        #### ERROR if number of reads are not the same (corrupt files)
            read_count +=1

            ## -- FORWARD READS -- ##

            # Read sequence and quality (encoded and decoded)
            read_fw = fastq_fw[1]
            qual_str_fw = fastq_fw[3]
            qual_score_fw = tf.quality_score(qual_str_fw, phred)

            ## Stats for original reads ##
            # Average length of all kept reads
            read_len_sum += len(read_fw)
            # Average quality of all kept reads
            avg_qual = sum(qual_score_fw)/len(qual_score_fw)
            read_qual_sum += avg_qual          

            ## STEP 0: Drop read if quality can't be determined
            if qual_score_fw == 'unknown':
                # Keep track of dropped read
                dropped_reads += 1
                # Reset
                line_count, fastq_fw, fastq_rev = 0, [], []
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline() 
                # Continue to next read without printing
                continue

            ## STEP 1: Remove leading and trailing bases, given user input
            read_fw, qual_str_fw, qual_score_fw = tf.removal_of_bases(read_fw, qual_str_fw, qual_score_fw, LEADING, TRAILING)

            ## STEP 2: Remove leading and trailing bases, based on quality 
            read_fw, qual_str_fw, qual_score_fw, trimmed = tf.sliding_window_pop(read_fw, qual_str_fw, qual_score_fw, WIN_SIZE, AVG_QUALITY, BASE_QUALITY)
            if (trimmed == True) and (trimmed_reads != 'all'):
                trimmed_reads += 1
            ## ------------------- ##
            
            ## -- REVERSE READS -- ##
            # Read sequence and quality (encoded and decoded)
            read_rev = fastq_rev[1]
            qual_str_rev = fastq_rev[3]
            qual_score_rev = tf.quality_score(qual_str_rev, phred)

            ## Stats for original reads ##
            # Average length of all kept reads
            read_len_sum += len(read_rev)
            # Average quality of all kept reads
            avg_qual = sum(qual_score_rev)/len(qual_score_rev)
            read_qual_sum += avg_qual            
            
            ## STEP 0: Drop read if quality can't be determined
            if qual_score_rev == 'unknown':
                # Keep track of dropped read
                dropped_reads += 1
                # Reset
                line_count, fastq_fw, fastq_rev = 0, [], []
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline()
                # Continue to next read without printing                
                continue

            ## STEP 1: Remove leading and trailing bases, given user input      
            read_rev, qual_str_rev, qual_score_rev = tf.removal_of_bases(read_rev, qual_str_rev, qual_score_rev, LEADING, TRAILING)

            ## STEP 2: Remove leading and trailing bases, based on quality
            read_rev, qual_str_rev, qual_score_rev, trimmed = tf.sliding_window_pop(read_rev, qual_str_rev, qual_score_rev, WIN_SIZE, AVG_QUALITY, BASE_QUALITY)
            if (trimmed == True) and (trimmed_reads != 'all'):
                trimmed_reads += 1
            ## ------------------- ##

            ## -- FORWARD AND REVERSE, SIMULTANEOUSLY -- ##
            ## STEP 3: Drop reads that become too short after trimming
            if (len(read_fw) < MIN_LEN) or (len(read_rev) < MIN_LEN):  
                # Keep track of dropped read
                dropped_reads += 1           
                # Reset
                line_count, fastq_fw, fastq_rev = 0, [], []
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline() 
                # Continue to next read without printing
                continue

            ## STEP 4: Drop reads with low average quality
            # Average quality in the reads
            avg_qual_fw = sum(qual_score_fw)/len(qual_score_fw)
            avg_qual_rev = sum(qual_score_rev)/len(qual_score_rev)
            # Drop trimmed reads with low average quality
            if (avg_qual_fw < AVG_QUALITY) or (avg_qual_rev < AVG_QUALITY):
                # Keep track of dropped read
                dropped_reads += 1                    
                # Reset
                line_count, fastq_fw, fastq_rev = 0, [], []
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline()    
                # Continue to next read without printing
                continue

            ## STEP 5: Drop reads with too many N bases
            if (read_fw.count('N') > N_MAX) or (read_rev.count('N') > N_MAX):
                # Keep track of dropped read
                dropped_reads += 1                               
                # Reset
                line_count, fastq_fw, fastq_rev = 0, [], []
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline()    
                # Continue to next read without printing
                continue 

            ## STEP 6: Print trimmed reads onto outfile, 
            if len(fastq_fw[0]) == 0 or len(fastq_rev[0]) == 0:
                print("Input files do not contain equal number of read. Output cannot be trused.")       ######### Added this, stops if number of reads differ
                sys.exit(1)

            tf.print_read(fastq_fw[0], read_fw, qual_str_fw, out_fw)
            tf.print_read(fastq_rev[0], read_rev, qual_str_rev, out_rev)
            
            ## Stats for trimmed reads
            # Average length of all kept reads
            trimmed_read_len_sum += len(read_fw) + len(read_rev)
            # Average quality of all kept reads
            trimmed_read_qual_sum += avg_qual_fw + avg_qual_rev
            
            ## ------------------- ##

            # Reset
            line_count, fastq_fw, fastq_rev = 0, [], []

        # Read next lines
        line_fw = file_fw.readline()
        line_rev = file_rev.readline()

       
    # Close
    file_fw.close()
    file_rev.close()
    out_fw.close()
    out_rev.close()

    ## ----------------------------- ##

    ## ---------- LOG FILE --------- ##
    log = open(base_fw + '.log', 'w')

    print('This is the log file for the trimming of', in_fwFile, 'and', in_revFile, file=log)
    print('\n===============\nSETTINGS\n===============', file=log)
    print('File 1:', in_fwFile, file=log)               # File 1
    print('File 2:', in_revFile, file=log)              # File 2
    print('Base quality: ', BASE_QUALITY, file=log)     # Single base quality
    print('Average quality:', AVG_QUALITY, file=log)    # Average quality
    print('Lead trim:', LEADING, file=log)              # Lead trim
    print('Trail trim:', TRAILING, file=log)            # Trail trim
    print('Window size:', WIN_SIZE, file=log)           # Window size
    print('Maximum unknown bases:', N_MAX, file=log)    # Maximum number of Ns
    print('Min lenght:', MIN_LEN, file=log)             # Min lenght after trim
    if (USER_PHRED != '') and (phred != USER_PHRED):    # Phred encoding type
        print("Phred encoding was set to {}, but {} was used.".format(USER_PHRED, phred), file=log)
    else:
        print('Phred: ' + phred, file=log)            

    print('\n===============\nSTATS\n===============', file=log)
    print('*** Before trimming ***', file=log)
    print('Total number of read pairs', read_count, file=log)
    print('Average length of kept reads:', round(read_len_sum/(read_count*2),2), file=log)
    print('Average quality of kept reads:', round(read_qual_sum/(read_count*2),2), file=log)
    print('\n*** After trimming ***', file=log)
    print('Read pairs dropped due to low quality/short length:', dropped_reads, file=log)
    print('Read pairs kept:', read_count-dropped_reads, file=log)
    print('Reads trimmed:', trimmed_reads, file=log)
    print('Average length of kept reads:', round(trimmed_read_len_sum/((read_count-dropped_reads)*2),2), file=log)
    print('Average quality of kept reads:', round(trimmed_read_qual_sum/((read_count-dropped_reads)*2),2), file=log)

    log.close()

    ## ------------------------- ##

    ## -------- STDOUT --------- ##
    # If trimming was successful
    print('Congratulations! Your trimming was successful. \nYou can find your results in the', base_fw + '_trimmed and', base_rev + '_trimmed files and some additional info in the the', base_fw + '.log file. \nPleasure working with you!')
    ## ------------------------- ## 