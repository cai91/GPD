def getGeneDensity(f):
    sizes=[]
    genes=[]
    p=''
    G=0
    with open(f) as inFile:
        for line in inFile:
            if line.strip()=='##FASTA': # Reached end of annotation
                break
            else:
                if line.strip()=='##gff-version 3':
                    continue
                if '##sequence-region' in line: # Read sizes of contigs
                    sizes.append(int(line.split()[-1]))
                else:
                    if line.split()[0]!=p:
                        genes.append(G)
                        p=line.split()[0]
                        G=0
                    G+=1
    genes.append(G)
    genes=genes[1:]
    
    return (np.array(genes)/(np.array(sizes)/1000.0))