def getSignature_hash(seq):
    
    k=5
    seq=seq.upper()
    
    all_bases=[]
    bases=['A','G','C','T']
    
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                for b4 in bases:
                    for b5 in bases:
                        all_bases.append(b1+b2+b3+b4+b5)
                          
    map_f={}
    for kmer in all_bases:
        map_f[kmer]=0

    # Count frequency of kmers in sequence      
    for i in range(0,len(seq)-(k-1)):
        try:
            map_f[seq[i:i+k]]+=1
        except:
            pass

    # Get normalised vector
    X=np.array([v for v in map_f.values()])
    X=X/sum(X)

    return X