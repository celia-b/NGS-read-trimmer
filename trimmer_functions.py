'''
Fucntions
*
'''
import argparse


def run_arg_parser():
    """An argparser function returning values from the command line."""
    parser = argparse.ArgumentParser(description = 'THE NGS READ TRIMMER')  # creating the argument parser

    parser.add_argument('file_name', metavar = 'File name 1',               # specifing argument (a positional argument)
                        help = "name of the fastq file")

    parser.add_argument('file_name2', metavar = 'File name 2',              # nargs='?' makes it an optinal positional argument 
                         nargs ='?' , default='',
                         help = "second file (optional)")

    parser.add_argument('-bq', '--base_qual', default='3', metavar='',
                        help="Quality threshold for single bases ")

    parser.add_argument('-aq', '--avg_qual', default='3', metavar='',
                        help="Average quality threshold for read ")

    parser.add_argument('-ls', '--lead_trim', default='0', metavar='',
                        help="X bases to trim from 3' strand")

    parser.add_argument('-ts', '--trail_trim', default='0', metavar='',
                        help="X bases to trim from 5' strand")

    parser.add_argument('-w', '--window_size', default='4', metavar='',     # optional arguments (metavar='' makes help output cleaner)
                        help='output to file bla bla')

    parser.add_argument('-t', '--threshold', default='15', metavar='',
                        help = 'quality score threshold')

    parser.add_argument('-ml', '--min_len', default='36', metavar='',
                        help="minimum length for trimmed reads")

    parser.add_argument('-ph', '--phred', default='', metavar='',
                        help="phred encoding, user input")

    return parser.parse_args() # getting the arguments from the parser


def quality_base_pop(read, qual_str, qual_score, BASE_QUALITY):
    """Removing 5 and 3 prime bases based on user input or default value"""
    leadingPopped = 0

    try:
        while qual_score[0] <= BASE_QUALITY:
            qual_score.pop(0)                                         # Removing low quality leading bases from translated quality string
            leadingPopped += 1                                        # Keep track of number of removed characters
    except IndexError:                                                # In case the entire read has low quality or all bases were removed in previous step 
        pass

    trailingPopped = 0
    try:
        while qual_score[-1] <= BASE_QUALITY:
            qual_score.pop(-1)                                        # Removing low quality leading bases from translated quality string
            trailingPopped += 1                                       # Keep track of number of removed characters
    except IndexError:                                                # In case the entire read has low quality or all bases were removed in previous step
        pass

    # Trim the read and the encoded quality string accordingly
    if trailingPopped != 0:
        read = read[leadingPopped:-trailingPopped] 
        qual_str = qual_str[leadingPopped:-trailingPopped]
    else: 
        read = read[leadingPopped:]
        qual_str = qual_str[leadingPopped:]

    return read, qual_str, qual_score



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
            return "phred33"

        elif not is_phred33 and is_phred64:
            phred_determined = True
            return "phred64"

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
    if phred == 'phred64':
        for i in range(len(quality_str)):
            quality_scores.append(int(ph64[quality_str[i]]))
        return quality_scores

    #If Phred33
    elif phred == 'phred33':
        for i in range(len(quality_str)):
            quality_scores.append(int(ph33[quality_str[i]]))
        return quality_scores

    # If string contains unknows characters None is returned
