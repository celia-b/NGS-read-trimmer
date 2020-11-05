def phred_control(fastq_file):
	"""Determining if a file is phred 33 or phred 64 encoded"""

	infile = open(fastq_file, 'r')
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

		line = infile.readline()[:-1]

	infile.close()

print(phred_control('test_small'))