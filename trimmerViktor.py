from trimmer_functions import run_arg_parser, phred_control, quality_score, removal_of_bases, quality_base_pop

args = run_arg_parser()

print('Settings\n===============')  
print('File 1: ' + args.file_name)          # File 1        (not used)
print('File 2: ' + args.file_name2)         # File 2        (not used)
print('Base quality: ' + args.base_qual)    # Single base quality
print('Average quality: ' + args.avg_qual)  # Average quality
print('Lead trim: ' + args.lead_trim)       # Lead trim
print('Trail trim: ' + args.trail_trim)     # Trail trim
print('Window size: ' + args.window_size)   # Window size   (not used)
print('Threshold: ' + args.threshold)       # Threshold     (not used)
print('Min lenght: ' + args.min_len)        # Min lenght after trim


LEADING = int(args.lead_trim)               # 3' bases to be removed
TRAILING = int(args.trail_trim)             # 5' bases to be removed
BASE_QUALITY = int(args.base_qual)          # Quality threshold for single bases 
AVG_QUALITY = int(args.avg_qual)            # Average quality threshold for read
MIN_LEN = int(args.min_len)                 # Minimum length for trimmed read

in_fwFile = args.file_name#'test_L2_1_pf.fastq'       ################### rename to file 1 and file 2?
in_rwFile = args.file_name2#'test_L2_2_pf.fastq'

#Singel read
############################ 
if in_rwFile == '':
    print("bara en fi fil")
    phred = phred_control(in_fwFile)  
    file_fw = open(in_fwFile, 'r')
    out_fw = open('trimmed_'+ in_fwFile, 'w')

    line_fw = file_fw.readline()
    line_count = 0
    fastq_fw = []

    while line_fw != '':
        line_count += 1
        fastq_fw.append(line_fw.strip())

        if line_count == 4:
            read_fw = fastq_fw[1]
            qual_str_fw = fastq_fw[3]
            qual_score_fw = quality_score(qual_str_fw, phred)
            # Remove leading and trailing bases, given user input
            read_fw, qual_str_fw, qual_score_fw = removal_of_bases(read_fw, qual_str_fw, qual_score_fw, LEADING, TRAILING)
            # Remove leading and trailing bases, based on quality 
            read_fw, qual_str_fw, qual_score_fw = quality_base_pop(read_fw, qual_str_fw, qual_score_fw, BASE_QUALITY)

            if (len(read_fw) < MIN_LEN):                                    
                line_count = 0
                fastq_fw = []
                line_fw = file_fw.readline()
                # Continue to next read without printing
                continue

            avgQualityForward = sum(qual_score_fw)/len(qual_score_fw)
            if (avgQualityForward < AVG_QUALITY):
                line_count = 0
                fastq_fw = []
                line_fw = file_fw.readline()
                # Continue to next read without printing
                continue

            print(fastq_fw[0], file=out_fw)
            print(read_fw, file=out_fw)
            print('+', file=out_fw)
            print(qual_str_fw, file=out_fw)

            # Reset
            line_count = 0
            fastq_fw = []

        line_fw = file_fw.readline()

    file_fw.close()
    out_fw.close()
    

#Pair end read
###############
else:
    #phred control, compare with file 2 if file 2 is included
    phred = phred_control(in_fwFile)  

    # opening files
    file_fw = open(in_fwFile, 'r')
    file_rev = open(in_rwFile, 'r')
    out_fw = open('trimmed_'+ in_fwFile, 'w')
    out_rev = open('trimmed_'+ in_rwFile, 'w')

    line_fw = file_fw.readline()
    line_rev = file_rev.readline() 

    line_count = 0

    fastq_fw = []
    fastq_rev = []

    while (line_fw or line_rev) != '':
        line_count += 1
        fastq_fw.append(line_fw.strip())
        fastq_rev.append(line_rev.strip())

        if line_count == 4:
        
            ## FORWARD READS ##
            read_fw = fastq_fw[1]
            qual_str_fw = fastq_fw[3]
            qual_score_fw = quality_score(qual_str_fw, phred)

            ## STEP 1: Remove leading and trailing bases, given user input
            read_fw, qual_str_fw, qual_score_fw = removal_of_bases(read_fw, qual_str_fw, qual_score_fw, LEADING, TRAILING)

            ## STEP 2: Remove leading and trailing bases, based on quality 
            read_fw, qual_str_fw, qual_score_fw = quality_base_pop(read_fw, qual_str_fw, qual_score_fw, BASE_QUALITY)

            ## REVERSE READS ##
            read_rev = fastq_rev[1]
            qual_str_rev = fastq_rev[3]
            qual_score_rev = quality_score(qual_str_rev, phred)

            ## STEP 1: Remove leading and trailing bases, given user input      
            read_rev, qual_str_rev, qual_score_rev = removal_of_bases(read_rev, qual_str_rev, qual_score_rev, LEADING, TRAILING)
          
            ## STEP 2: Remove leading and trailing bases, based on quality
            read_rev, qual_str_rev, qual_score_rev = quality_base_pop(read_rev, qual_str_rev, qual_score_rev, BASE_QUALITY)

            ## STEP 3: Drop reads that become too short after trimming
            if (len(read_fw) < MIN_LEN) or (len(read_rev) < MIN_LEN):                                    

                # Reset
                line_count = 0
                fastq_fw = []
                fastq_rev = []
                
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline() 
                # Continue to next read without printing
                continue

            ## STEP 4: Drop reads with low average quality
            # Average quality in the forward read

            avgQualityForward = sum(qual_score_fw)/len(qual_score_fw)

            # Average quality in the reverse read
        
            avgQualityReverse = sum(qual_score_rev)/len(qual_score_rev)

            # Drop trimmed reads with low average quality
            if (avgQualityForward < AVG_QUALITY) or (avgQualityReverse < AVG_QUALITY):
                # Reset
                line_count = 0
                fastq_fw = []
                fastq_rev = []

                #reset_lines(line_count, fastq_fw, fastq_rev) ################################ ?????????????
                
                # Read new lines
                line_fw = file_fw.readline()
                line_rev = file_rev.readline()    

                # Continue to next read without printing
                continue


            # STEP 5: Print trimmed reads onto outfile 

            print(fastq_fw[0], file=out_fw)
            print(read_fw, file=out_fw)
            print('+', file=out_fw)
            print(qual_str_fw, file=out_fw)

            print(fastq_rev[0], file=out_rev)
            print(read_rev, file=out_rev)
            print('+', file=out_rev)
            print(qual_str_rev, file=out_rev)


            # Reset
            line_count = 0
            fastq_fw = []
            fastq_rev = []

        line_fw = file_fw.readline()
        line_rev = file_rev.readline()

    ## -------------- ##


    file_fw.close()
    file_rev.close()
    out_fw.close()
    out_rev.close()




