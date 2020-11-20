import trimmer_functions as tf
import gzip

args = tf.run_arg_parser()

## Print settings
print('===============\nSETTINGS\n===============')  
print('File 1: ' + args.file_name)          # File 1
print('File 2: ' + args.file_name2)         # File 2
print('Base quality: ' + args.base_qual)    # Single base quality
print('Average quality: ' + args.avg_qual)  # Average quality
print('Lead trim: ' + args.lead_trim)       # Lead trim
print('Trail trim: ' + args.trail_trim)     # Trail trim
print('Window size: ' + args.window_size)   # Window size   (not used)
print('Threshold: ' + args.threshold)       # Threshold     (not used)
print('Min lenght: ' + args.min_len)        # Min lenght after trim
print('Phred: ' + args.phred)

## Take parameters and filenames from arguments
LEADING = int(args.lead_trim)               # 3' bases to be removed
TRAILING = int(args.trail_trim)             # 5' bases to be removed
BASE_QUALITY = int(args.base_qual)          # Quality threshold for single bases 
AVG_QUALITY = int(args.avg_qual)            # Average quality threshold for read
MIN_LEN = int(args.min_len)                 # Minimum length for trimmed read
N_MAX = 3
WIN_SIZE = int(args.window_size)
USER_PHRED = 'phred33'

in_fwFile = args.file_name
in_revFile = args.file_name2

log = open('BRISCOE_log.txt', 'w')

##################
# Singe end read #
################## 
if in_revFile == '':
    # Determine Phred encoding type
    phred = tf.phred_control(in_fwFile, USER_PHRED) 
    
    # Open files, checking if file is compressed or not 
    if in_fwFile.endswith('.gz'):
        file_fw = gzip.open(in_fwFile, 'rt')
    else:    
        file_fw = open(in_fwFile, 'r')
    
    # Output file, add '.gz' file ending if uncompressed file was read
    if in_fwFile.endswith('.gz'):
        out_fw = gzip.open('trimmed_' + in_fwFile, 'wt')
    else:    
        out_fw = gzip.open('trimmed_' + in_fwFile + '.gz', 'wt')

    # Read first line and initialize variables
    line_fw = file_fw.readline()
    
    line_count = 0
    fastq_fw = []

    dropped_reads = 0

    # Iterate through lines
    while line_fw != '':
        line_count += 1
        # Make a 4 element list with the 4 lines of each read
        fastq_fw.append(line_fw.strip())

        if line_count == 4:
            # Read sequence
            read_fw = fastq_fw[1]

            # Read quality: encoded (str) and decoded (score)
            qual_str_fw = fastq_fw[3]

            qual_score_fw = tf.quality_score(qual_str_fw, phred)

            # If quality score cannt be calculate, don't print to output file
            if qual_score_fw == 'unknown':
                continue
            
            ## STEP 1: Remove leading and trailing bases, given user input
            read_fw, qual_str_fw, qual_score_fw = tf.removal_of_bases(read_fw, qual_str_fw, qual_score_fw, LEADING, TRAILING)
            
            ## STEP 2: Remove leading and trailing bases, based on quality
            read_fw, qual_str_fw, qual_score_fw = tf.sliding_window_pop(read_fw, qual_str_fw, qual_score_fw, WIN_SIZE, AVG_QUALITY)

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
            avgQualityForward = sum(qual_score_fw)/len(qual_score_fw)
            if (avgQualityForward < AVG_QUALITY):
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
            print(fastq_fw[0], file=out_fw)
            print(read_fw, file=out_fw)
            print('+', file=out_fw)
            print(qual_str_fw, file=out_fw)

            # Reset
            line_count, fastq_fw = 0, []

        # Read next line
        line_fw = file_fw.readline()

    # Close files
    file_fw.close()
    out_fw.close()



####################
# Paired end reads #
####################
else:

    if tf.phred_control(in_fwFile, USER_PHRED) != tf.phred_control(in_revFile, USER_PHRED):                     ###### Controlling the files use the same phred 
        print("ERROR")

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

    # Output file, add '.gz' file ending if uncompressed file was read
    if in_fwFile.endswith('.gz'):
        out_fw = gzip.open('trimmed_' + in_fwFile, 'wt')
    else:    
        out_fw = gzip.open('trimmed_' + in_fwFile + '.gz', 'wt')
        
    if in_revFile.endswith('.gz'):
        out_rev = gzip.open('trimmed_' + in_revFile, 'w')
    else:    
        out_rev = gzip.open('trimmed_' + in_revFile + '.gz', 'wt')
        

    # Read first line of each file and initialize variables
    line_fw = file_fw.readline()
    line_rev = file_rev.readline() 

    line_count, fastq_fw, fastq_rev = 0, [], []

    dropped_reads = 0

    # Iterate through lines
    while (line_fw or line_rev) != '':
        line_count += 1
        
        # Make a 4 element list with the 4 lines of each read, per file
        fastq_fw.append(line_fw.strip())
        fastq_rev.append(line_rev.strip())

        if line_count == 4:
        
            ## FORWARD READS ##
            # Read sequence
            read_fw = fastq_fw[1]
            
            # Read quality: encoded (str) and decoded (score)
            qual_str_fw = fastq_fw[3]
            qual_score_fw = tf.quality_score(qual_str_fw, phred)

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
            read_fw, qual_str_fw, qual_score_fw = tf.sliding_window_pop(read_fw, qual_str_fw, qual_score_fw, WIN_SIZE, AVG_QUALITY)

            
            ## REVERSE READS ##
            # Read sequence
            read_rev = fastq_rev[1]
            
            # Read quality: encoded (str) and decoded (score)
            qual_str_rev = fastq_rev[3]

            qual_score_rev = tf.quality_score(qual_str_rev, phred)
            
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
            read_rev, qual_str_rev, qual_score_rev = tf.sliding_window_pop(read_rev, qual_str_rev, qual_score_rev, WIN_SIZE, AVG_QUALITY)
            
            ## STEP 3 (COMMON): Drop reads that become too short after trimming
            if (len(read_fw) < MIN_LEN) or (len(read_rev) < MIN_LEN):  
                # Keep track of dropped read
                dropped_reads += 1
                # Print to log file
                print('REMOVED READ!', file=log)
                print('This read is', len(read_fw), 'bases long', file=log)
                print(fastq_fw[0], file=log)
                print(read_fw, file=log)
                print('+', file=log)
                print(qual_str_fw, file=log)
                
                print('Its pair is', len(read_rev), 'bases long', file=log)
                print(fastq_rev[0], file=log)
                print(read_rev,file=log)
                print('+', file=log)
                print(qual_str_rev, file=log) 
                print('', file=log)                 
                # Reset
                line_count, fastq_fw, fastq_rev = 0, [], []
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline() 
                # Continue to next read without printing
                continue

            ## STEP 4 (COMMON): Drop reads with low average quality
            # Average quality in the forward read
            avgQualityForward = sum(qual_score_fw)/len(qual_score_fw)
            # Average quality in the reverse read
            avgQualityReverse = sum(qual_score_rev)/len(qual_score_rev)
            # Drop trimmed reads with low average quality
            if (avgQualityForward < AVG_QUALITY) or (avgQualityReverse < AVG_QUALITY):
                # Keep track of dropped read
                dropped_reads += 1
                # Print to log file
                print('REMOVED READ!', file=log)
                print('This read is', len(read_fw), 'bases long', file=log)
                print(fastq_fw[0], file=log)
                print(read_fw, file=log)
                print('+', file=log)
                print(qual_str_fw, file=log)
                
                print('Its pair is', len(read_rev), 'bases long', file=log)
                print(fastq_rev[0], file=log)
                print(read_rev,file=log)
                print('+', file=log)
                print(qual_str_rev, file=log) 
                print('', file=log)                                 
                # Reset
                line_count, fastq_fw, fastq_rev = 0, [], []
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline()    
                # Continue to next read without printing
                continue

            # STEP 5: Drop reads with too many N bases
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

            ## STEP 5 (COMMON): Print trimmed reads onto outfile 
            print(fastq_fw[0], file=out_fw)
            print(read_fw, file=out_fw)
            print('+', file=out_fw)
            print(qual_str_fw, file=out_fw)

            print(fastq_rev[0], file=out_rev)
            print(read_rev, file=out_rev)
            print('+', file=out_rev)
            print(qual_str_rev, file=out_rev)


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
    log.close()


## Print stats
print('')
print('===============\nSTATS\n===============')
print('Read pairs dropped due to low quality/short length:', dropped_reads)

if phred != USER_PHRED:
    print("Phred encoding was set to {}, but {} was used.".format(USER_PHRED, phred))                   ######  If input phred is not the same as determined phred.