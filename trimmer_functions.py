'''
Functions
'''
import argparse
import gzip
import sys
import os

def run_arg_parser():
    """An argparser function returning values from the command line."""
    parser = argparse.ArgumentParser(description = 'THE NGS READ TRIMMER')  # creating the argument parser

    parser.add_argument('FILE1', metavar = 'File name 1',               # specifing argument (a positional argument)
                        help = "Your forward fastq file.")

    parser.add_argument('FILE2', metavar = 'File name 2',              # nargs='?' makes it an optional positional argument 
                         nargs ='?' , default='',
                         help = "Your reverse fastq file (optional, only for paired-end reads).")
    
    parser.add_argument('-PH',ﬁ '--PHRED', default='', metavar='',
                        help = " Phred encoding type (phred 33 or phred 64). The encoding type will automatically be determined by the program, and if user input does not match true type, a warning will be printed onto the log file.")

    parser.add_argument('-L', '--LEADING', default='0', metavar='',
                        help = "Number of bases to be trimmed from 3' end of all reads, regardless of quality. Default is 0.")

    parser.add_argument('-T', '--TRAILING', default='0', metavar='',
                        help = "Number of bases to be trimmed from all 5' end of all reads, regardless of quality. Default is 0.")

    parser.add_argument('-W', '--WINDOWSIZE', default='4', metavar='',   
                        help = 'Window size for sliding window trimming approach. Default is 4. If WINDOWSIZE is 1, then the single base appoach is used and BASEQUALITY is used as quality threshold, instead of AVGQUALITY.')

    parser.add_argument('-AQ', '--AVGQUALITY', default='15', metavar='',
                        help = "Average quality threshold for sliding window trimming approach. Default is 15.")

    parser.add_argument('-BQ', '--BASEQUALITY', default='3', metavar='',
                        help = "Quality threshold for single base trimming approach. Default is 3.")

    parser.add_argument('-ML', '--MINLEN', default='36', metavar='',
                        help = "Minimum allowed length of reads. Default is 36.")

    parser.add_argument('-N', '--MAXN', default='15', metavar='',
                        help = 'Maximum number of unknown bases allowed in a read. Default is 3.')

    return parser.parse_args() # getting the arguments from the parser


def sliding_window_pop(read, qual_str, qual_score, WIN_SIZE, AVG_QUALITY, BASE_QUALITY):
    """Removing 5 and 3 prime bases with sliding window approach"""
    trimmed = False
    ## LEADING TRIM 
    # Read first window
    window = qual_score[:WIN_SIZE]
    # Slide  
    if WIN_SIZE > 1:      
        while (len(window) != 0) and (sum(window)/len(window) < AVG_QUALITY):
            trimmed = True
            # Remove 3' end base
            qual_score = qual_score[1:]
            qual_str = qual_str[1:]
            read = read[1:]
            # Read new window
            window = qual_score[:WIN_SIZE]
    elif WIN_SIZE == 1:
        while (len(window) != 0) and (sum(window)/len(window) < BASE_QUALITY):
            trimmed = True
            # Remove 3' end base
            qual_score = qual_score[1:]
            qual_str = qual_str[1:]
            read = read[1:]
            # Read new window
            window = qual_score[:WIN_SIZE]   

    ## TRAILING TRIM
    # Read first window
    window = qual_score[-WIN_SIZE:]
    # Slide
    if WIN_SIZE > 1:
        while (len(window) != 0) and (sum(window)/len(window) < AVG_QUALITY):
            trimmed = True
            # Remove 5' end base
            qual_score = qual_score[:-1]
            qual_str = qual_str[:-1]
            read = read[:-1]
            # Read new window
            window = qual_score[-WIN_SIZE:]
    elif WIN_SIZE == 1:
        while (len(window) != 0) and (sum(window)/len(window) < BASE_QUALITY):
            trimmed = True
            # Remove 5' end base
            qual_score = qual_score[:-1]
            qual_str = qual_str[:-1]
            read = read[:-1]
            # Read new window
            window = qual_score[-WIN_SIZE:]        

    return read, qual_str, qual_score, trimmed


def removal_of_bases(DNA_str,quality_str, quality_score, LEADING, TRAILING): 
    """Remove leading and trailing bases, given user input"""
    if TRAILING != 0:
        quality_score = quality_score[LEADING:-TRAILING]
        quality_str = quality_str[LEADING:-TRAILING]
        DNA_str = DNA_str[LEADING:-TRAILING]
    else: 
        quality_score = quality_score[LEADING:]
        quality_str = quality_str[LEADING:]
        DNA_str = DNA_str[LEADING:]

    return DNA_str, quality_str, quality_score


def phred_control(fastqFile, user_phred):
    """Determining if a file is phred 33 or phred 64 encoded"""

    if fastqFile.endswith('.gz'):
        infile = gzip.open(fastqFile, 'rt')
    else:    
        infile = open(fastqFile, 'r')

    line = infile.readline()[:-1]

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
            return "33"

        elif not is_phred33 and is_phred64:
            phred_determined = True
            return "64"

        # In case phred can't be determined, use the users input. 
        if line == '': 
            phred_determined = True
            return user_phred

        line = infile.readline()[:-1]

    infile.close()


def quality_score(quality_str, phred):

    # Controlling if line contains unknown character 
    set_ph64 = set("@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh")
    set_ph33 = set("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ")
    qual_set = set(quality_str)

    if len(qual_set.difference(set_ph33))>0 and len(qual_set.difference(set_ph64))>0:
        return 'unknown'
    
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
    if phred == '64':
        for i in range(len(quality_str)):
            quality_scores.append(int(ph64[quality_str[i]]))
        return quality_scores

    #If Phred33
    elif phred == '33':
        for i in range(len(quality_str)):
            quality_scores.append(int(ph33[quality_str[i]]))
        return quality_scores

    # If string contains unknows characters None is returned


def print_read(ID, seq, qual_str, file):
    print(ID, file=file)
    print(seq, file=file)
    print('+', file=file)
    print(qual_str, file=file)

def controling_output_file(in_file):
    if in_file.endswith('.gz'):
        base = in_file.split('.')[0]
        file_name = base + '_trimmed.txt.gz'
        if  os.path.exists(os.getcwd() + '/' + file_name):
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
        base = in_file.split('.')[0]
        file_name = base + '_trimmed.fa'
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