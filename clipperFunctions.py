''' -----------------------------------------
    These are the functions needed for the magicClipper NGS read trimmer.
    -----------------------------------------

    Please have this file along with magicClipper.py in your desired directory for
    correct functioning.
    
    -----------------------------------------
    Authors: Celia Burgos Sequeros (s202423) & Viktor Törnblom (s200116)
    Last update: November 2020
    -----------------------------------------
'''

## Required modules
import argparse
import gzip
import sys
import os

def run_arg_parser():
    """
    This is an argparser function that returns arguments given in the command line.
    """
    
    parser = argparse.ArgumentParser(description = '---- THE magicClipper NGS READ TRIMMER ---- \
    When given one or two .fastq files, this program performs a complete user-guided \
    and quality-based trimming of the reads contained in it. \
    If you like, you can adjust the settings listed below. \
    Otherwise, their defaults will be used.')

    parser.add_argument('FILE1', metavar = 'File name 1',
                        help = "Your forward fastq file.")

    parser.add_argument('FILE2', metavar = 'File name 2',              # nargs='?' makes it an optional positional argument 
                         nargs ='?' , default='',
                         help = "Your reverse fastq file (optional, only for paired end mode).")
    
    parser.add_argument('-PH', '--PHRED', default='', metavar='',
                        help = " Phred encoding type (phred 33 or phred 64). \
                        The encoding type will automatically be determined by the program, \
                        and if user input does not match true type, a warning will be printed \
                        onto the log file.")

    parser.add_argument('-L', '--LEADING', default='0', metavar='',
                        help = "Number of bases to be trimmed from 3' end of all reads, \
                        regardless of quality. Default is 0.")

    parser.add_argument('-T', '--TRAILING', default='0', metavar='',
                        help = "Number of bases to be trimmed from 5' end of all reads, \
                        regardless of quality. Default is 0.")

    parser.add_argument('-W', '--WINDOWSIZE', default='4', metavar='',   
                        help = 'Window size for sliding window trimming approach. Default is 4. \
                        If WINDOWSIZE is 1, then the single base appoach is used and BASEQUALITY \
                        is taken to be the quality threshold, instead of AVGQUALITY.')

    parser.add_argument('-AQ', '--AVGQUALITY', default='15', metavar='',
                        help = "Average quality threshold for sliding window trimming approach. \
                        Default is 15.")

    parser.add_argument('-BQ', '--BASEQUALITY', default='3', metavar='',
                        help = "Quality threshold for single base trimming approach. Default is 3.")

    parser.add_argument('-ML', '--MINLEN', default='36', metavar='',
                        help = "Minimum allowed length of reads. Default is 36.")

    parser.add_argument('-N', '--MAXN', default='15', metavar='',
                        help = 'Maximum number of unknown bases allowed in a read. Default is 15.')

    return parser.parse_args()


def quality_trim(read, qual_str, qual_score, WIN_SIZE, AVG_QUALITY, BASE_QUALITY):
    """
    This function removes 5' and 3' based on their quality.
    If WIN_SIZE > 1, it takes the sliding window approach, with AVG_QUALITY as threshold.
    If WIN_SIZE = 1, it takes the single base approach, with BASE_QUALITY as threshold.
    """
    
    trimmed = False                         # With this flag we keep track of the trimmed reads
    
    ## LEADING TRIM 
    window = qual_score[:WIN_SIZE]          # Read first window 
    if WIN_SIZE > 1:                        # Slide     
        while (len(window) != 0) and (sum(window)/len(window) < AVG_QUALITY):
            trimmed = True
            
            qual_score = qual_score[1:]     # Remove 3' end base
            qual_str = qual_str[1:]
            read = read[1:]
            
            window = qual_score[:WIN_SIZE]  # Read new window
    
    elif WIN_SIZE == 1:                     # Slide
        while (len(window) != 0) and (sum(window)/len(window) < BASE_QUALITY):
            trimmed = True
            
            qual_score = qual_score[1:]     # Remove 3' end base
            qual_str = qual_str[1:]
            read = read[1:]
            
            window = qual_score[:WIN_SIZE]  # Read new window   

    ## TRAILING TRIM
    window = qual_score[-WIN_SIZE:]         # Read first window
    if WIN_SIZE > 1:                        # Slide
        while (len(window) != 0) and (sum(window)/len(window) < AVG_QUALITY):
            trimmed = True
            
            qual_score = qual_score[:-1]    # Remove 5' end base
            qual_str = qual_str[:-1]
            read = read[:-1]
            
            window = qual_score[-WIN_SIZE:] # Read new window
    
    elif WIN_SIZE == 1:                     # Slide
        while (len(window) != 0) and (sum(window)/len(window) < BASE_QUALITY):
            trimmed = True
            
            qual_score = qual_score[:-1]    # Remove 5' end base
            qual_str = qual_str[:-1]
            read = read[:-1]
            
            window = qual_score[-WIN_SIZE:] # Read new window      

    return read, qual_str, qual_score, trimmed


def global_trim(DNA_str,quality_str, quality_score, LEADING, TRAILING): 
    """
    This function removes leading and trailing bases in all reads, if user says so.
    """
    if TRAILING != 0:
        quality_score = quality_score[LEADING:-TRAILING]
        quality_str = quality_str[LEADING:-TRAILING]
        DNA_str = DNA_str[LEADING:-TRAILING]
    
    else: 
        quality_score = quality_score[LEADING:]
        quality_str = quality_str[LEADING:]
        DNA_str = DNA_str[LEADING:]

    return DNA_str, quality_str, quality_score


def phred_autodetect(input_file, USER_PHRED):
    """
    This function detects if a file is encoded using phred33 or phred64
    """

    if input_file.endswith('.gz'):              # Open file
        infile = gzip.open(input_file, 'rt')
    else:    
        infile = open(input_file, 'r')               

    # Phred sets
    phred64_set = set("@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh")
    phred33_set = set("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ")

    quality_string = ''                         # Initialize variables
    line_count = 0
    is_phred33 = False                          
    is_phred64 = False
    phred_determined = False

    line = infile.readline()[:-1]               # Read line by line, until phred type is found
    while phred_determined == False:
        line_count += 1

        if line_count == 4:                     # At this point, we are looking at a quality string
            quality_string = line
            is_phred33 = set(quality_string).issubset(phred33_set)
            is_phred64 = set(quality_string).issubset(phred64_set)
            line_count = 0

            if is_phred33 and not is_phred64:
                phred_determined = True
                return "33"

            elif not is_phred33 and is_phred64:
                phred_determined = True
                return "64"
        
        line = infile.readline().strip()

    infile.close()

    # In case phred can't be determined, use the users input. 
    if not phred_determined: 
        # If user did not specify phred type 
        if USER_PHRED == '':
            print('ERROR: We cannot autodetect the phred encoding type of your file(s). Please specify it in the input.')
            sys.exit(1)
        phred_determined = True
        return USER_PHRED  


def quality_score(quality_str, phred):

    # Controlling if line contains unknown character 
    phred64_set = set("@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh")
    phred33_set = set("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ")
    qual_set = set(quality_str)

    phred64_dict = {'@': 0,'A': 1,'B': 2,'C': 3, 'D': 4,'E': 5,'F': 6,'G': 7,'H': 8,'I': 9,
        'J': 10,'K': 11,'L': 12,'M': 13,'N': 14,'O': 15,'P': 16,'Q': 17,'R': 18,
        'S': 19,'T': 20,'U': 21,'V': 22,'W': 23,'X': 24,'Y': 25,'Z': 26,'[': 27,
        '\\': 28,']': 29,'^': 30,'_': 31,'`': 32,'a': 33,'b': 34,'c': 35,'d': 36,
        'e': 37, 'f': 38,'g': 39,'h': 40,'i': 41}

    phred33_dict = {'!': 0,'"': 1,'#': 2,'$': 3,'%': 4,'&': 5,"'": 6,'(': 7,')': 8,
        '*': 9,'+': 10,',': 11,'-': 12,'.': 13,'/': 14,'0': 15,'1': 16,
        '2': 17,'3': 18,'4': 19,'5': 20,'6': 21,'7': 22,'8': 23,'9': 24,
        ':': 25,';': 26,'<': 27,'=': 28,'>': 29,'?': 30,'@': 31,'A': 32,
        'B': 33,'C': 34,'D': 35,'E': 36,'F': 37,'G': 38,'H': 39,'I': 40,
        'J': 41}

    quality_scores = []  
    
    if phred == '64':
        for i in range(len(quality_str)):
            try:
                quality_scores.append(int(phred64_dict[quality_str[i]]))
            except KeyError:                        # If unknown character is found, return 'unknown'. This read will be removed in the main body of the program.
                quality_scores='unknown'
                break
        return quality_scores

    #If Phred33
    elif phred == '33':
        for i in range(len(quality_str)):
            try:
                quality_scores.append(int(phred33_dict[quality_str[i]]))
            except KeyError:
                quality_scores='unknown'
                break
        return quality_scores  
  

def print_read(ID, seq, qual_str, file):
    print(ID, file=file)
    print(seq, file=file)
    print('+', file=file)
    print(qual_str, file=file)


def controling_output_file(input_file):
    if input_file.endswith('.gz'):
        base = input_file.split('.')[0]
        file_name = base + '_trimmed.fastq.gz'
        if os.path.exists(os.getcwd() + '/' + file_name):
            answer = None
            while  answer not in ['y','n']:
                answer = input("{} will be overwritten. Do you want to continue? y/n: ".format(file_name))
                if answer == 'y':
                    return(gzip.open(file_name, 'wt'))
                elif answer == 'n':
                    print('Exiting program')
                    sys.exit(1)
                else:
                    print('Invalid input')
        else:
            return(gzip.open(file_name, 'wt'))
    else:
        base = input_file.split('.')[0]
        file_name = base + '_trimmed.fastq'
        if  os.path.exists(os.getcwd() + '/' + file_name):
            answer = None
            while  answer not in ['y','n']:
                answer = input("{} will be overwritten. Do you want to continue? y/n: ".format(file_name))
                if answer == 'y':
                    return(open(file_name, 'w'))
                elif answer == 'n':
                    print('Exiting program')
                    sys.exit(1)
                else:
                    print('Invalid input')
        else:
            return(open(file_name, 'w'))