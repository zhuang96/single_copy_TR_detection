# NOT TESTED; skeleton code suggested by Marzie

info=[]
line=0
i=0 

#format every 60 bases
def printer():
	new_chr21.write("\n")
	line=0

def reading_info(info_input):  
    first_ind=int(info_input[0])
    last_ind=int(info_input[1])
    pattern_size=int(info_input[2])
    copy_num=float(info_input[3])
    fasta_header=info_input[5]
    arr_seq=info_input[6]
    
    final_len= pattern_size*(copy_num - (trunc(copy_num-1.8)+1))
    return final_len

def trunc(x):
    if (x<1):
        return 0
    else:
        return np.floor(x)


#open TRDB file
with open('/Users/Zoe/Desktop/SPRING18/EXTRA_CIR/TANDEMREPEATS/ref_set_chr21.txt') as t, open("/Users/Zoe/Desktop/SPRING18/EXTRA_CIR/TANDEMREPEATS/chr21.fa") as f:
	    fa_chr21=open("fa_chr21.fa","w+")
	    fa_chr21.write(">mod_chr21\n") #chrX
	    idx_chr21=open("idx_chr21.fa","w+")
	    idx_chr21.write("FASTA\tstart\tend\n")
	    next(f)
	    next(t) 
	    tline=next(t)
	    for fline in f:
	        for c in fline:
	        	info = tline.split("\t")
	        	#while right of TR 
	        	while (i < start):
	        		fa_chr21.write(str(c))
	        		i+=1
	        		line+=1
	        		c+=1
	        		if (line==60):
	        			printer()
	        	#get the length of the truncated TR
	        	trun_len=reading_info(info)
	        	#write FASTA header, start index
	        	idx_chr21.write(str(info[5])+"\t"+str(i)+"\t")
	        	for x in trun_len:
	        		#write start of new index to other file too
	        		fa_chr21.write(str(c))
	        		# idx_chr21.write(str(c))
	        		i+=1
	        		c+=1
	        		line+=1
	        		if (line==60):
	        			printer()
	        	#write end index
	        	idx_chr21.write(str(i)+"\n")
	        	#if end of ref file
	        	if next(t) not t:
	        		#while not end of .fa file
	        		while (c in f):
	        			fa_chr21.write(str(c))
	        			idx_chr21.close()
	        			fa_chr21.close()
	        			break
	        	#update start 
	        	tline=next(t)

