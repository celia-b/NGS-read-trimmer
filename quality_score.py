def quality_score(quality_str):
    
    # the set of characters in the quality line from file 
    quality_characters = set(quality_str)
    phred64 = set("@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh")
    phred33 = set("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ")
    
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
    
    # If not possible to distinguish between phred33 and phred64
    if quality_characters.issubset(phred64) and quality_characters.issubset(phred33):
        print('BOTH! ERROR!')
        
    elif quality_characters.issubset(phred64):
        for i in range(len(qual1)):
            quality_scores.append(int(ph64[qual1[i]]))
    
    elif quality_characters.issubset(phred33):
        for i in range(len(qual1)):
             quality_scores.append(int(ph33[qual1[i]]))
   
    # If string contains unknows acsii characters
    else:
        print("Error")
    
    return quality_scores

qual1 = "#1=DDFFFHHFHHJJIJJJJJJJJGJIJJJJIJIIJIIJJJJJJIJJJJJJIJJJJJJJJJJJJHHHHFFFFFEEAEEDDDDDDDDDDDDDDDCDDDDDDD"
print(quality_score(qual1))
