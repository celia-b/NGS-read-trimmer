'''
Fucntions
*
*
*

'''
import argparse


def run_arg_parser():
	"""An argparser function returning values from the command line."""
	parser = argparse.ArgumentParser(description = 'THE NGS READ TRIMMER') 	# creating the argument parser

	parser.add_argument('file_name', metavar = 'File name 1',				# specifing argument (a positional argument)
						help = "name of the fastq file")

	parser.add_argument('file_name2', metavar = 'File name 2',				# nargs='?' makes it an optinal positional argument 
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

	parser.add_argument('-w', '--window_size', default='4', metavar='', 	# optional arguments (metavar='' makes help output cleaner)
						help='output to file bla bla')

	parser.add_argument('-t', '--threshold', default='15', metavar='',
						help = 'quality score threshold')

	parser.add_argument('-ml', '--min_len', default='36', metavar='',
						help="minimum length for trimmed reads")

	return parser.parse_args() # getting the arguments from the parser




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




def quality_score(quality_str, phred):

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