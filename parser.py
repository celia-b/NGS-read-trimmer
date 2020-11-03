import argparse

parser = argparse.ArgumentParser(description = 'THE NGS READ TRIMMER') 	# creating the argument parser

parser.add_argument('file_name', metavar = 'File name 1',				# specifing argument (a positional argument)
					help = "name of the fastq file")

parser.add_argument('file_name2', metavar = 'File name 2',				# nargs='?' makes it an optinal positional argument 
					 nargs ='?' , default='',
                     help = "second file (optional)")

parser.add_argument('-w', '--window_size', default='4', metavar='', 	# optional arguments (metavar='' makes help output cleaner)
					help='output to file bla bla')

parser.add_argument('-t', '--threshold', default='15', metavar='',
					help = 'quality score threshold')

args = parser.parse_args() # getting the arguments from the parser

print('File name: ' + args.file_name)
print('File name: ' + args.file_name2)
print('Window size: ' + args.window_size)
print('Threshold: ' + args.threshold)

###############################################################################################

'''
infile1 = open(args.file_name, 'r')
infile2 = open(args.file_name2, 'r')

line1 = infile1.readline()
line2 = infile2.readline() 

line_count = 0
fastq1 = []
fastq2 = []

while (line1 or line2) != '':
    line_count += 1
    fastq1.append(line1[:-1])
    fastq2.append(line2[:-1])
    
    if line_count == 4:
        print(fastq1)
        print()
        print(fastq2)
        print()
        print()
        line_count = 0
        fastq1 = []
        fastq2 = []
        
    line1 = infile1.readline()
    line2 = infile2.readline()
'''