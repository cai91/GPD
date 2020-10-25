def getHprot(f):
    sizes=[]
    genes=[]
    hprots=[]
    p=''
    G=0
    H=0
    with open(f) as inFile:
        for line in inFile:
            if line.strip()=='##FASTA': # Reached end of annotation
                break
            else:
                if line.strip()=='##gff-version 3' or '##sequence-region' in line:
                    continue
                else:
                    if line.split()[0]!=p:
                        genes.append(G)
                        hprots.append(H)
                        p=line.split()[0]
                        
                        G=0
                        H=0
                        
                    G+=1
                    # Detect hypothetical proteins
                    for t in line.split('\t')[8].split(';'):
                        if t.split('=')[0]=='product':
                            prod=t.split('=')[1]
                            if prod=='hypothetical protein':
                                H+=1
    genes.append(G)
    hprots.append(H)
    
    genes=genes[1:]
    hprots=hprots[1:]
    
    
    return np.array(hprots)/np.array(genes)