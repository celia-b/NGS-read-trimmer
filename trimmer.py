from trimmer_functions import run_arg_parser, phred_control, quality_score

args = run_arg_parser()

'''
print('Settings\n===============')	
print('File 1: ' + args.file_name)			# File 1 		(not used)
print('File 2: ' + args.file_name2)			# File 2 		(not used)
print('Base quality: ' + args.base_qual)	# Single base quality
print('Average quality: ' + args.avg_qual)	# Average quality
print('Lead trim: ' + args.lead_trim)		# Lead trim
print('Trail trim: ' + args.trail_trim)		# Trail trim
print('Window size: ' + args.window_size)	# Window size 	(not used)
print('Threshold: ' + args.threshold)		# Threshold 	(not used)
print('Min lenght: ' + args.min_len)		# Min lenght after trim
'''

LEADING = int(args.lead_trim)				# 3' bases to be removed
TRAILING = int(args.trail_trim)				# 5' bases to be removed
BASE_QUALITY = int(args.base_qual)			# Quality threshold for single bases 
AVG_QUALITY = int(args.avg_qual) 			# Average quality threshold for read
MIN_LEN = int(args.min_len)					# Minimum length for trimmed read

in_fwFile = 'test_L2_1_pf.fastq' 	# remove later and replace with argparser input
in_rwFile = 'test_L2_2_pf.fastq'

# opening files
fileForward = open(in_fwFile, 'r')
fileReverse = open(in_rwFile, 'r')
outForward = open('trimmed_'+ in_fwFile, 'w')
outReverse = open('trimmed_'+ in_rwFile, 'w')

#phred control, compare with file 2 if file 2 is included
phred = phred_control(fileForward)	


lineForward = fileForward.readline()
lineReverse = fileReverse.readline() 

line_count = 0

fastqForward = []
fastqReverse = []





def removal_of_bases(DNA_str,quality_str, quality_score, LEADNING, TRAILING):    
    if TRAILING != 0:
        quality_score = quality_score[LEADING:-TRAILING]
        quality_str = quality_str[LEADING:-TRAILING]
        read = read[LEADING:-TRAILING]
    else: 
        quality_score = quality_score[LEADING:]
        quality_str = quality_str[LEADING:]
        read = read[LEADING:]

    return read, quality_str, quality_score

readForward, qualityStringForward, qualityScoreForward = removal_of_bases(readForward, qualityStringForward, qualityScoreForward, LEADING, TRAILING)

readFw, qualStrFw, qualScoreFw = removal_of_bases(readFw, qualStrFw, qualScoreFw, LEADING, TRAILING)



while (lineForward or lineReverse) != '':									########## We also need to be able to deal with one file
    line_count += 1

    fastqForward.append(lineForward.strip())
    fastqReverse.append(lineReverse.strip())

    if line_count == 4:
        ## FORWARD READS ##
        readForward = fastqForward[1]
        qualityStringForward = fastqForward[3]
        qualityScoreForward = quality_score(qualityStringForward, phred)

        ## STEP 1: Remove leading and trailing bases, given user input
        if TRAILING != 0:
            qualityScoreForward = qualityScoreForward[LEADING:-TRAILING]
            qualityStringForward = qualityStringForward[LEADING:-TRAILING]
            readForward = readForward[LEADING:-TRAILING]
        else: 
            qualityScoreForward = qualityScoreForward[LEADING:]
            qualityStringForward = qualityStringForward[LEADING:]
            readForward = readForward[LEADING:]      

        ## STEP 2: Remove leading and trailing bases, based on quality       ########### Rewrite to moving window?
        # Remove leading bases
        leadingPopped = 0

        try:
            while qualityScoreForward[0] <= BASE_QUALITY:
                qualityScoreForward.pop(0)                                # Removing low quality leading bases from translated quality string
                leadingPopped += 1                                        # Keep track of number of removed characters
        except IndexError:                                                # In case the entire read has low quality or all bases were removed in previous step 
            pass

        # Remove trailing bases
        trailingPopped = 0
        try:
            while qualityScoreForward[-1] <= BASE_QUALITY:
                qualityScoreForward.pop(-1)                               # Removing low quality leading bases from translated quality string
                trailingPopped += 1                                       # Keep track of number of removed characters
        except IndexError:                                                # In case the entire read has low quality or all bases were removed in previous step
            pass

        # Trim the read and the encoded quality string accordingly
        if trailingPopped != 0:
            readForward = readForward[leadingPopped:-trailingPopped] 
            qualityStringForward = qualityStringForward[leadingPopped:-trailingPopped]
        else: 
            readForward = readForward[leadingPopped:]
            qualityStringForward = qualityStringForward[leadingPopped:]



        ## REVERSE READS ##
        readReverse = fastqReverse[1]
        qualityStringReverse = fastqReverse[3]
        qualityScoreReverse = quality_score(qualityStringReverse, phred)

        ## STEP 1: Remove leading and trailing bases, given user input
        if TRAILING != 0:
            qualityScoreReverse = qualityScoreReverse[LEADING:-TRAILING]
            qualityStringReverse = qualityStringReverse[LEADING:-TRAILING]
            readReverse = readReverse[LEADING:-TRAILING]
        else: 
            qualityScoreReverse = qualityScoreReverse[LEADING:]
            qualityStringReverse = qualityStringReverse[LEADING:]
            readReverse = readReverse[LEADING:]


        ## STEP 2: Remove leading and trailing bases, based on quality
        # Remove leading bases
        leadingPopped = 0

        try:
            while qualityScoreReverse[0] <= BASE_QUALITY:
                qualityScoreReverse.pop(0)                                # Removing low quality leading bases from translated quality string
                leadingPopped += 1                                        # Keep track of number of removed characters
        except IndexError:                                                # In case the entire read has low quality or all bases were removed in previous step
            pass                                                          

        # Remove trailing bases
        trailingPopped = 0
        try: 
            while qualityScoreReverse[-1] <= BASE_QUALITY:
                qualityScoreReverse.pop(-1)                               # Removing low quality leading bases from translated quality string
                trailingPopped += 1                                       # Keep track of number of removed characters
        except IndexError:                                                # In case the entire read has low quality or all bases were removed in previous step
            pass

        # Trim the read and the encoded quality string accordingly
        if trailingPopped != 0:
            readReverse = readReverse[leadingPopped:-trailingPopped] 
            qualityStringReverse = qualityStringReverse[leadingPopped:-trailingPopped]
        else: 
            readReverse = readReverse[leadingPopped:]
            qualityStringReverse = qualityStringReverse[leadingPopped:]


        ## STEP 3: Drop reads that become too short after trimming
        if (len(readForward) < MIN_LEN) or (len(readReverse) < MIN_LEN):                                    
            # Reset
            line_count = 0
            fastqForward = []
            fastqReverse = [] 

            # Read new lines
            lineForward = fileForward.readline()
            lineReverse = fileReverse.readline() 

            # Continue to next read without printing
            continue


        ## STEP 4: Drop reads with low average quality
        # Average quality in the forward read
        sumQualityForward = 0
        for elem in qualityScoreForward:									######### sum(qualityScoreSorward)
            sumQualityForward += elem
        avgQualityForward = sumQualityForward/len(qualityScoreForward)

        # Average quality in the reverse read
        sumQualityReverse = 0
        for elem in qualityScoreReverse:
            sumQualityReverse += elem
        avgQualityReverse = sumQualityReverse/len(qualityScoreReverse)

        # Drop trimmed reads with low average quality
        if (avgQualityForward < AVG_QUALITY) or (avgQualityReverse < AVG_QUALITY):
            # Reset
            line_count = 0
            fastqForward = []
            fastqReverse = []

            # Read new lines
            lineForward = fileForward.readline()
            lineReverse = fileReverse.readline()    

            # Continue to next read without printing
            continue


        # STEP 5: Print trimmed reads onto outfile						
        print(fastqForward[0], file=outForward)
        print(readForward, file=outForward)
        print('+', file=outForward)
        print(qualityStringForward, file=outForward)

        print(fastqReverse[0], file=outReverse)
        print(readReverse, file=outReverse)
        print('+', file=outReverse)
        print(qualityStringReverse, file=outReverse)


        # Reset
        line_count = 0
        fastqForward = []
        fastqReverse = []

    lineForward = fileForward.readline()
    lineReverse = fileReverse.readline()

## -------------- ##


fileForward.close()
fileReverse.close()
outForward.close()
outReverse.close()




