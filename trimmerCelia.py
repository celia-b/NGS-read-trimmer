## -- My files -- ##
#fileForward = open('testfile_1.fastq', 'r')
#fileReverse = open('testfile_2.fastq', 'r')

#outForward = open('outtest_1.fastq', 'w')
#outReverse = open('outtest_2.fastq', 'w')


fileForward = open('BRISCOE_0069_BD18RUACXX_L2_1_pf.fastq', 'r')
fileReverse = open('BRISCOE_0069_BD18RUACXX_L2_2_pf.fastq', 'r')

outForward = open('trimmed_BRISCOE_0069_BD18RUACXX_L2_1_pf.fastq', 'w')
outReverse = open('trimmed_BRISCOE_0069_BD18RUACXX_L2_2_pf.fastq', 'w')

## -------------- ##


## -- My default parameters -- ##
# Leading and trailing bases to be removed
LEADING = 0
TRAILING = 0

# Quality threshold for single bases 
BASE_QUALITY = 3

# Average quality threshold for read
AVG_QUALITY = 15

# Minimum length for trimmed read
MIN_LEN = 36
## --------------------------- ##


## -- My functions -- ##
def phred_control(fastqFile):
	"""Determining if a file is phred 33 or phred 64 encoded"""

	line = fastqFile.readline()[:-1]

	phred64 = set("@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh")
	phred33 = set("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ")

	quality_string = ''
	line_count = 0
	is_phred33 = False
	is_phred64 = False
	phred_determined = False

	while phred_determined == False:
		line_count += 1

		if line_count == 4:
			quality_string = line
			is_phred33 = set(quality_string).issubset(phred33)
			is_phred64 = set(quality_string).issubset(phred64)
			line_count = 0

		if is_phred33 and not is_phred64:
			phred_determined = True
			return "phred33"

		elif not is_phred33 and is_phred64:
			phred_determined = True
			return "phred64"

		line = fastqFile.readline()[:-1]

def quality_score(quality_str):
    
    ph64 = {'@': 0,'A': 1,'B': 2,'C': 3, 'D': 4,'E': 5,'F': 6,'G': 7,'H': 8,'I': 9,
        'J': 10,'K': 11,'L': 12,'M': 13,'N': 14,'O': 15,'P': 16,'Q': 17,'R': 18,
        'S': 19,'T': 20,'U': 21,'V': 22,'W': 23,'X': 24,'Y': 25,'Z': 26,'[': 27,
        '\\': 28,']': 29,'^': 30,'_': 31,'`': 32,'a': 33,'b': 34,'c': 35,'d': 36,
        'e': 37, 'f': 38,'g': 39,'h': 40,'i': 41}

    ph33 = {'!': 0,'"': 1,'#': 2,'$': 3,'%': 4,'&': 5,"'": 6,'(': 7,')': 8,
        '*': 9,'+': 10,',': 11,'-': 12,'.': 13,'/': 14,'0': 15,'1': 16,
        '2': 17,'3': 18,'4': 19,'5': 20,'6': 21,'7': 22,'8': 23,'9': 24,
        ':': 25,';': 26,'<': 27,'=': 28,'>': 29,'?': 30,'@': 31,'A': 32,
        'B': 33,'C': 34,'D': 35,'E': 36,'F': 37,'G': 38,'H': 39,'I': 40,
        'J': 41}
    
    quality_scores = []

    # If Phred64   
    if phred == 'phred64':
        for i in range(len(quality_str)):
            quality_scores.append(int(ph64[quality_str[i]]))
    
    #If Phred33
    elif phred == 'phred33':
        for i in range(len(quality_str)):
            quality_scores.append(int(ph33[quality_str[i]]))
   
    # If string contains unknows characters
    else:
        print("There are unknown characters in the quality string of your read. This read will be ommited from the analysis.")
    
    return quality_scores
## ------------------ ##


## -- Check PHRED encoding -- ## (To make it more robust we should probably check both files..?)
phred = phred_control(fileForward)

# Need to close and re-open file bc otherwise parser will get confused
fileForward.close()
fileReverse.close()
fileForward = open('BRISCOE_0069_BD18RUACXX_L2_1_pf.fastq', 'r')
fileReverse = open('BRISCOE_0069_BD18RUACXX_L2_2_pf.fastq', 'r')
## -------------------------- ##



## -- Do the trimming -- ##
# We will be working with fastq files. Each read in a fastqc file looks like this:

# HEADER LINE:      @HWI-ST1294:69:D18RUACXX:2:1101:1160:2048_1:N:0:
# SEQUENCE LINE:    NGCCAGACCCTATCAGAAAGAAAAGAAAAGAGAAGAGAAGAGAAGAGGAAAGGAAAGAGAAAGAGAAGGAAAGAAAGAAAGAAAGAGAAGGAAAGAAAGAA
# SEPARATOR LINE:   +
# QUALITY LINE:     #1=DDFFFHHFHHJJIJJJJJJJJGJIJJJJIJIIJIIJJJJJJIJJJJJJIJJJJJJJJJJJJHHHHFFFFFEEAEEDDDDDDDDDDDDDDDCDDDDDDD


lineForward = fileForward.readline()
lineReverse = fileReverse.readline() 

line_count = 0

fastqForward = []
fastqReverse = []

while (lineForward or lineReverse) != '':
    line_count += 1
    
    fastqForward.append(lineForward.strip())
    fastqReverse.append(lineReverse.strip())

    if line_count == 4:
        ## FORWARD READS ##
        readForward = fastqForward[1]
        qualityStringForward = fastqForward[3]
        qualityScoreForward = quality_score(qualityStringForward)

        ## STEP 1: Remove leading and trailing bases, given user input
        if TRAILING != 0:
            qualityScoreForward = qualityScoreForward[LEADING:-TRAILING]
            qualityStringForward = qualityStringForward[LEADING:-TRAILING]
            readForward = readForward[LEADING:-TRAILING]
        else: 
            qualityScoreForward = qualityScoreForward[LEADING:]
            qualityStringForward = qualityStringForward[LEADING:]
            readForward = readForward[LEADING:]      
        
        ## STEP 2: Remove leading and trailing bases, based on quality
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
        qualityScoreReverse = quality_score(qualityStringReverse)

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
        for elem in qualityScoreForward:
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
