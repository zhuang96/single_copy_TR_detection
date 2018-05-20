import numpy as np

#global variables
ref_list=[]

fasta_list=[]
chrstart_list=[]
chrend_list=[]
check_tr=[]


#truncate and add to ref_list the FIRST_INDEX, NEW_COPY, LAST_INDEX, FASTA_HEADER
def reading_info(info_input):
    
    first_ind=int(info_input[0])
    last_ind=int(info_input[1])
    pattern_size=int(info_input[2])
    copy_num=float(info_input[3])
    fasta_header=info_input[5]
    arr_seq=info_input[6]
    
    final_len= pattern_size*(copy_num - (trunc(copy_num-1.8)+1))

    single_copy=arr_seq[:int(round(final_len))]

    ref_list.append([first_ind,single_copy,last_ind,fasta_header]) 

def trunc(x):
    if (x<1):
        return 0
    else:
        return np.floor(x)
    
def create_new_seq(ref_string):
    new_ref=""
    #tracker of starting index
    start_ind=0
    
    # fasta_list.append([ref_list[0][-1]])
    # chrstart_list.append([start_ind]) #wrong 

    #go through oredered_dictionary to append flanking left seq and TR single copy sequence
    for i in range(len(ref_list)):
        # if ([ref_list[i][-1]] not in fasta_list):
        #     #print ref_list[i][-1]
        #     fasta_list.append([ref_list[0][-1]])
        #     chrstart_list.append([start_ind-1])
        #     chrend_list.append(len(new_ref)-1)

        # print(ref_string[start_ind:cur_ind])
        # print("start_ind is " + str(start_ind) + " cur_ind is " + str(cur_ind))
        # if(i==2):
        # 	break;


        fasta_list.append([ref_list[i][-1]])
        #start of TR 
        cur_ind=ref_list[i][0]
        #append sequence w/o TR
        new_ref+=(ref_string[start_ind:cur_ind])
        #start of TR is len(new_ref) b/c +1 -1 prior start index
        chrstart_list.append([len(new_ref)-10])
        #append shortened TR copy
        new_ref+=ref_list[i][1]
        #end of TR is len(new_ref)-1
        chrend_list.append([len(new_ref)-1+10])
        #double check TR indexes
        check_tr.append(new_ref[chrstart_list[i][0]:chrend_list[i][0]+1])
        #start again at end index
        start_ind=ref_list[i][2]
        
    #append last chunk of string after last detected TR
    new_ref+=(ref_string[start_ind:])

    # chrend_list.append([len(new_ref)-1])

    return new_ref
    
#def main():
#open TRDB file
with open('/Users/Zoe/Desktop/SPRING18/EXTRA_CIR/TANDEMREPEATS/ref_set_chr21.txt') as f:
    #skip first line, then do
    next(f) 
    for line in f:
        reading_info(line.split("\t"))
#open ref seq file 
with open("/Users/Zoe/Desktop/SPRING18/EXTRA_CIR/TANDEMREPEATS/chr21.fa") as dat:
    #skip first line
    next(dat)
    #since the file is based 0, must add a white space char
    ref_nowhitespace=" "
    for line in dat:
        #remove all whitespace
        ref_nowhitespace += line.strip() 
        
#make whole sequence UPPERCASE
ref_nowhitespace=ref_nowhitespace.upper()


create_modified_chr21=open("modified_chr21.txt","w+")
create_modified_chr21.write(ref_nowhitespace)
create_modified_chr21.close()

#testing ; add one because also want the last index
# print ref_nowhitespace[5010564:5010609+1]
#    print ref_nowhitespace[46699842:46699983+1]
#    print ref_nowhitespace[46506725:46506819+1]

new_ref_string=create_new_seq(ref_nowhitespace)

################################TESTING####################################

# for i in range(len(ref_list)):
# 	print str(ref_list[i])
# print new_ref_string[5010564+2454:5010564+2454+35]
# print ref_list[1]
# print ref_list[-8]
# print len(ref_list[-8][1])


################################CREATING FILES####################################
#double checking 
# dc=open("chr21_feed.fa","w+")
# for i in range(len (fasta_list)):
# 	dc.write(str(fasta_list[i][0])+"\t"+str(chrstart_list[i][0]-1)+"\t"+str(chrend_list[i][0]-1)+"\n")
#     # dc.write(str(fasta_list[i][0])+"\t"+str(chrstart_list[i][0])+"\t"+str(chrend_list[i][0])+"\t"+str(check_tr[i])+"\t"+str(ref_list[i][1]+"\n"))
# dc.close() 

# checkthis=open("checkthis.txt","w+")
# sindex=45787976
# for i in range(4036,4041):
# 	checkthis.write("prior: "+ str(new_ref_string[sindex:chrstart_list[i][0]]+"\n"+"compare to: "+str(ref_list[i])))
# 	sindex=chrstart_list[i][0]
# checkthis.close()
# create genome file
g=open("chr21_new.bed","w+")
# g.write("FASTA_HEADER\tSTART\tEND\n")
for i in range(len (fasta_list)):
    g.write(str(fasta_list[i][0])+"\t"+str(chrstart_list[i][0])+"\t"+str(chrend_list[i][0])+"\n")
g.close() 

# print len (fasta_list)

# create wraparound of new string
# wraparound= open("chr21_new_wraparound.fa","w+")
# wraparound.write(">chr21 new_wraparound")
# wraparound.write(new_ref_string)
# wraparound.close() 

# #create wraparound of old string
# wraparound= open("chr21_wraparound.txt","w+")
# wraparound.write(ref_nowhitespace)
# wraparound.close() 

