# MagicClipper - An NGS Read Trimmer
The advent of Next Generation Sequencing (NGS) technologies have transformed how biological research is being performed and today almost all biological fields use the technology for cutting edge discoveries. From sequencing a human genome to studying whole bacterial communities and their interplay with the environment, the possibilities that NGS represents are unprecedented. However, these experiments produce massive amounts of data that often including a lot of noise and the program therefore we must to be able to clean it up. 

The tool that is developed in this project takes care of a crucial initial step in NGS data processing: removing low-quality reads and read-ends to later minimize noise in the alignment. The program:
  1. Reads and writes both compressed and uncompressed fastq files according to user input.
  2. Can understand both Phred+33 and Phred+64 encoding according to user input, or autodetects the scale if none is given.
  3. Trims nucleotides from the 3' and 5' ends of each read given user input.
  4. Trims each read from the 3' and 5' ends based on quality, calculated as mean of moving window.
  5. Filters out reads with a mean quality or length lower than specified after trimming, as well as those that have a more than a specified number of unknown bases.
  6. Keeps track of trimmed/removed reads and outputs a summary at the end of the analysis.

The program is written in Python 3.7. and can be run on any NGS file in FASTQ format. A more detailed explaination of the program, how to run it and what output to expect can be found in the file MagicClipper Report in this repository. 
