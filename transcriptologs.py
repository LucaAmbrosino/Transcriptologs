''' 
MIT License

Copyright (c) 2016 LucaAmbrosino

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

'''

import time
g=time.time()
import shutil
import os
import argparse

parser = argparse.ArgumentParser(description='Transcriptologs: Process the results from two reciprocal tBLASTx all-vs-all to predict orthologs.')

parser.add_argument("-i1", "--input1", type=str, action='store', dest='input1file', help="Input the 1st tBLASTx results file") 
parser.add_argument("-i2", "--input2", type=str, action='store', dest='input2file', help="Input the 2nd tBLASTx results file") 
parser.add_argument("-o", "--output", action='store', dest='outputfile', help="Directs the output to a name of your choice")
args = parser.parse_args()




os.mkdir(".transcriptologsTEMP")

f1=open(args.input1file)
f2=open(".transcriptologsTEMP\\output1_part1.txt","w")

n=0

for line in f1:
    if line.startswith("Query="):
        
        newline=line.split()
        query= newline[1]
        #print (query)
    elif line.startswith("Length="):
        query_length= line[7:-1]
        
    elif line.startswith(">"):
        f2.write("Query"+"\t"+"Subject"+"\t"+"Q. length"+"\t"+"S. length"+"\t"+"Q. start"+"\t"+"Q. end"+"\t"+"S. start"+"\t"+"S. end"
         +"\t"+"Identities"+"\t"+"Positives"+"\t"+"Gaps"+"\t"+"Frame"+"\t"+"Score"+"\t"+"E-value"+"\n")
        subjelements=line.split()
        subject= subjelements[1]
        line=f1.readline()
        while 1:
            if line.startswith("Length="):
                subject_length=line[7:-1]
                break
            line=f1.readline()

    elif line.startswith(" Score ="):
        #print(line)
                        
        elements1=line.split()
        score=elements1[2]
        evalue=elements1[7]
        line=f1.readline()
        elements2=line.split()
        identities=elements2[2]+" "+elements2[3][0:-1]
        positives=elements2[6]+elements2[7][0:-1]
        gaps=elements2[10]+elements2[11]
        line=f1.readline()
        elements3=line.split()
        frame=elements3[2]
        line=f1.readline()
        line=f1.readline()
        n=0
        
        while line.startswith("Query "):
            fragments1=line.split()
            if n==0:
                query_start=fragments1[1]
            query_end=fragments1[3]
            line=f1.readline()
            line=f1.readline()
            fragments2=line.split()
            if n==0:
                sbj_start=fragments2[1]
            sbj_end=fragments2[3]
            line=f1.readline()
            line=f1.readline()
            n=n+1
            
        f2.write(query+"\t"+subject+"\t"+query_length+"\t"+subject_length+"\t"+query_start+"\t"+
                query_end+"\t"+sbj_start+"\t"+sbj_end+"\t"+identities+"\t"+positives+"\t"+gaps+"\t"+frame+"\t"+score+"\t"+evalue+"\n")
                        
            
f2.write("Query")  
f1.close()
f2.close()

f1=open(".transcriptologsTEMP\\output1_part1.txt")
f2=open(".transcriptologsTEMP\\output1_part1_TEMP.txt","w")
f3=open(".transcriptologsTEMP\\output1_part2.txt","w")

templist=[]
templist2=[]
c=0
score=0
fragments=0
positives=0
identities=0
gaps=0
alignment_len=0
ref="vuoto"
for line in f1:
    if not line.startswith("Query"):
        elements=tuple(line.split())
        
        
        templist.append(elements)
    elif line.startswith("Query") and len(templist)!=0:
        f2.write("Query"+"\t"+"Subject"+"\t"+"Q. length"+"\t"+"S. length"+"\t"+"Q. start"+"\t"+"Q. end"+"\t"+"S. start"+"\t"+"S. end"
         +"\t"+"Identities"+"\t"+"Positives"+"\t"+"Gaps"+"\t"+"Frame"+"\t"+"Score"+"\t"+"E-value"+"\n")
        #templist.sort(key=lambda obj:float(obj[13]), reverse=True)
        templist2.append(templist[0])
        reference_frame=templist[0][12]
        reference_evalue=templist[0][14]
        
   
        for tupla in templist[1:]:
            if tupla[12][0]==reference_frame[0] and tupla[12][3]==reference_frame[3]:
                c=0


                for tupla2 in templist2:
                
                    #f3.write(str(tupla)+"\n"+str(tupla2)+"\n"+"\n")
                    if tupla2[12][0]=="+" and tupla2[12][3]=="+":
                        
                        if int(tupla[4])>=(int(tupla2[5])-20) or int(tupla[5])<=(int(tupla2[4])+20):
                            if int(tupla[6])>=(int(tupla2[7])-20) or int(tupla[7])<=(int(tupla2[6])+20):
                                c=c+1
                                pass
                    if tupla2[12][0]=="+" and tupla2[12][3]=="-":
                        if int(tupla[4])>=(int(tupla2[5])-20) or int(tupla[5])<=(int(tupla2[4])+20):
                            if int(tupla[6])<=(int(tupla2[7])+20) or int(tupla[7])>=(int(tupla2[6])-20):
                                c=c+1
                                pass
                    if tupla2[12][0]=="-" and tupla2[12][3]=="-":
                        if int(tupla[4])<=(int(tupla2[5])+20) or int(tupla[5])>=(int(tupla2[4])-20):
                            if int(tupla[6])<=(int(tupla2[7])+20) or int(tupla[7])>=(int(tupla2[6])-20):
                                c=c+1
                                pass
                    if tupla2[12][0]=="-" and tupla2[12][3]=="+":
                        if int(tupla[4])<=(int(tupla2[5])+20) or int(tupla[5])>=(int(tupla2[4])-20):
                            if int(tupla[6])>=(int(tupla2[7])-20) or int(tupla[7])<=(int(tupla2[6])+20):
                                c=c+1
                                pass
                    if c==len(templist2):
                        templist2.append(tupla)
                        c=0                       
                        
                
                    #f3.write(str(tupla2)+"\n")
                    
                    
        for tupla2 in templist2:
            elements_for_identities=tupla2[8].split("/")
            identities=identities+int(elements_for_identities[0])
            elements_for_positives=tupla2[10].split("/")
            positives=positives+int(elements_for_positives[0])
            new_elem=elements_for_positives[1].split("(")
            alignment_len=alignment_len+int(new_elem[0])
            elements_for_gaps=tupla2[11].split("/")
            gaps=gaps+int(elements_for_gaps[0])
            score=score+float(tupla2[13])
            fragments=fragments+1
            query_length=int(tupla2[2])
            subject_length=int(tupla2[3])
            
            for element in tupla2:
                
                f2.write(element+"\t")
            f2.write("\n")
        if templist[0][0]!=ref:
            f3.write("# query"+"\t"+"subject"+"\t"+"query_length"+"\t"+"subject_length"+"\t"+"aln_legth"+"\t"+"identities"+"\t"+"positives"+"\t"+"gaps"+"\t"+"frame"+"\t"+"score"+"\t"+"n._of_fragments"+"\t"+"E-Value"+"\n")
        ref=templist[0][0]
        #print (query_length/3, alignment_len, subject_length/3)
        #if int(alignment_len)<=(query_length/3) and alignment_len<=(subject_length/3):    
        f3.write(templist[0][0]+"\t"+templist[0][1]+"\t"+str(round(query_length/3))+"\t"+str(round(subject_length/3))+"\t"+str(alignment_len)+"\t"+str(identities)+"\t"+str(positives)+"\t"+str(gaps)+"\t"+tupla2[12]+"\t"+str(round(score))+"\t"+str(fragments)+"\t"+reference_evalue+"\n")
        #elif int(alignment_len)>(query_length/3):
        #    f3.write(templist[0][0]+"\t"+templist[0][1]+"\t"+str(round(score))+"\t"+str(fragments)+"\t"+str(round(query_length/3))+"\t"+str(round(subject_length/3))+"\t"+str(round(query_length/3))+"\t"+str(identities)+"\t"+str(float(identities/query_length/3))+"\t"+str(positives)+"\t"+str(float(positives/query_length/3))+"\t"+str(gaps)+"\n")
        #elif int(alignment_len)>(subject_length/3):
        #    f3.write(templist[0][0]+"\t"+templist[0][1]+"\t"+str(round(score))+"\t"+str(fragments)+"\t"+str(round(query_length/3))+"\t"+str(round(subject_length/3))+"\t"+str(round(subject_length/3))+"\t"+str(identities)+"\t"+str(float(identities/subject_length/3))+"\t"+str(positives)+"\t"+str(float(positives/subject_length/3))+"\t"+str(gaps)+"\n")
        score=0
        fragments=0
        identities=0
        gaps=0
        positives=0
        alignment_len=0
        templist2=[]
        templist=[]
f1.close()
f2.close()
f3.close()

f1=open(".transcriptologsTEMP\\output1_part2.txt")
#f2=open("C:\Documents and Settings\Luca\Desktop\\provaBEST_HIT3.txt", "w")
f2=open(".transcriptologsTEMP\\output1_part3.txt", "w")

templist=[]
templist2=[]

for line in f1:
    if not line.startswith("#"):
        templist.append(line)
        
    elif line.startswith("#") and len(templist)!=0:
        for item in templist:
            elements=item.split()
            score=float(elements[9])
                
                #f2.write(elements[0]+"\t"+elements[1]+"\t"+elements[10]+"\t"+str(round(totalscore))+"\n")
            templist2.append((round(score),elements[0],elements[1],elements[2],elements[3],elements[4],elements[5],elements[6],elements[7],elements[8],elements[10],elements[11]))
                
        a=sorted(templist2,reverse=True)
                
        reference=a[0]
                
        ref_score=float(reference[0])
                
        for element in templist2:
                    
            current_score=element[0]
            if current_score >= ref_score:
                    
                f2.write(element[1]+"\t"+element[2]+"\t"+str(element[3])+"\t"+str(element[4])+"\t"+str(element[5])+"\t"+str(element[6])+"\t"+str(element[7])+"\t"+str(element[8])+"\t"+str(element[9])+"\t"+str(element[0])+"\t"+str(element[10])+"\t"+str(element[11])+"\n")
                    
        templist=[]
        templist2=[]
f1.close()
f2.close()

f1=open(".transcriptologsTEMP\\output1_part3.txt")

f2=open(".transcriptologsTEMP\\output1_part3_TEMP.txt", "w")

ALL=f1.readlines()
ALL.sort()
for line in ALL:
    f2.write(line)
f2.close()
f1.close()
f2=open(".transcriptologsTEMP\\output1_part3_TEMP.txt")
f3=open(".transcriptologsTEMP\\output1_part4.txt", "w")

p1="vuoto"
p2="vuoto"
p3=100000
previousline="primalinea"
templist=[]
for line in f2:
    elements=tuple(line.split())
    n1=elements[0]
    
    n2=elements[1]
    n3=elements[2]
    if len(templist)==0:
        templist.append(elements)
    
    elif not n1==p1:
        
        templist.sort(key=lambda obj:int(obj[2]), reverse=True)
        reference=templist[0][2]

        
        for item in templist:
            
            if int(item[2])>=int(reference):
                f3.write(item[0]+"\t"+item[1]+"\t"+item[2]+"\t"+item[3]+"\t"+item[4]+"\t"+item[5]+"\t"+item[6]+"\t"+item[7]+"\t"+item[8]+"\t"+item[9]+"\t"+item[10]+"\t"+item[11]+"\n")
        templist=[]
        templist.append(elements)
    
    elif n1==p1 and n2==p2 and n3!=p3:
        templist.append(elements)
    elif n1==p1 and n2!=p2:
        templist.append(elements)
        
        
    p1=elements[0]
    p2=elements[1]
    p3=elements[2]
f1.close()
f2.close()
f3.close()

f1=open(".transcriptologsTEMP\\output1_part4.txt")

f2=open(".transcriptologsTEMP\\output1_part5.txt", "w")

ALL=[]

for line in f1:
    elements=line.split()
    newelements=elements[0]+"\t"+elements[1]+"\t"+elements[2]+"\t"+elements[3]+"\t"+elements[4]+"\t"+elements[5]+"\t"+elements[6]+"\t"+elements[7]+"\t"+elements[8]+"\t"+elements[9]+"\t"+elements[10]+"\t"+elements[11]+"\n"
    ALL.append(newelements)

ALL.sort()
for element in ALL:
    f2.write(element)
f1.close()
f2.close()

f1=open(args.input2file)
f2=open(".transcriptologsTEMP\\output2_part1.txt","w")

n=0

for line in f1:
    if line.startswith("Query="):
        
        newline=line.split()
        query= newline[1]
        #print (query)
    elif line.startswith("Length="):
        query_length= line[7:-1]
        
    elif line.startswith(">"):
        f2.write("Query"+"\t"+"Subject"+"\t"+"Q. length"+"\t"+"S. length"+"\t"+"Q. start"+"\t"+"Q. end"+"\t"+"S. start"+"\t"+"S. end"
         +"\t"+"Identities"+"\t"+"Positives"+"\t"+"Gaps"+"\t"+"Frame"+"\t"+"Score"+"\t"+"E-value"+"\n")
        subjelements=line.split()
        subject= subjelements[1]
        line=f1.readline()
        while 1:
            if line.startswith("Length="):
                subject_length=line[7:-1]
                break
            line=f1.readline()

    elif line.startswith(" Score ="):
        #print(line)
                        
        elements1=line.split()
        score=elements1[2]
        evalue=elements1[7]
        line=f1.readline()
        elements2=line.split()
        identities=elements2[2]+" "+elements2[3][0:-1]
        positives=elements2[6]+elements2[7][0:-1]
        gaps=elements2[10]+elements2[11]
        line=f1.readline()
        elements3=line.split()
        frame=elements3[2]
        line=f1.readline()
        line=f1.readline()
        n=0
        
        while line.startswith("Query "):
            fragments1=line.split()
            if n==0:
                query_start=fragments1[1]
            query_end=fragments1[3]
            line=f1.readline()
            line=f1.readline()
            fragments2=line.split()
            if n==0:
                sbj_start=fragments2[1]
            sbj_end=fragments2[3]
            line=f1.readline()
            line=f1.readline()
            n=n+1
            
        f2.write(query+"\t"+subject+"\t"+query_length+"\t"+subject_length+"\t"+query_start+"\t"+
                query_end+"\t"+sbj_start+"\t"+sbj_end+"\t"+identities+"\t"+positives+"\t"+gaps+"\t"+frame+"\t"+score+"\t"+evalue+"\n")
                        
            
f2.write("Query")  
f1.close()
f2.close()

f1=open(".transcriptologsTEMP\\output2_part1.txt")
f2=open(".transcriptologsTEMP\\output2_part1_TEMP.txt","w")
f3=open(".transcriptologsTEMP\\output2_part2.txt","w")

templist=[]
templist2=[]
c=0
score=0
fragments=0
positives=0
identities=0
gaps=0
alignment_len=0
ref="vuoto"
for line in f1:
    if not line.startswith("Query"):
        elements=tuple(line.split())
        
        
        templist.append(elements)
    elif line.startswith("Query") and len(templist)!=0:
        f2.write("Query"+"\t"+"Subject"+"\t"+"Q. length"+"\t"+"S. length"+"\t"+"Q. start"+"\t"+"Q. end"+"\t"+"S. start"+"\t"+"S. end"
         +"\t"+"Identities"+"\t"+"Positives"+"\t"+"Gaps"+"\t"+"Frame"+"\t"+"Score"+"\t"+"E-value"+"\n")
        #templist.sort(key=lambda obj:float(obj[13]), reverse=True)
        templist2.append(templist[0])
        reference_frame=templist[0][12]
        reference_evalue=templist[0][14]
        
   
        for tupla in templist[1:]:
            if tupla[12][0]==reference_frame[0] and tupla[12][3]==reference_frame[3]:
                c=0


                for tupla2 in templist2:
                
                    #f3.write(str(tupla)+"\n"+str(tupla2)+"\n"+"\n")
                    if tupla2[12][0]=="+" and tupla2[12][3]=="+":
                        
                        if int(tupla[4])>=(int(tupla2[5])-20) or int(tupla[5])<=(int(tupla2[4])+20):
                            if int(tupla[6])>=(int(tupla2[7])-20) or int(tupla[7])<=(int(tupla2[6])+20):
                                c=c+1
                                pass
                    if tupla2[12][0]=="+" and tupla2[12][3]=="-":
                        if int(tupla[4])>=(int(tupla2[5])-20) or int(tupla[5])<=(int(tupla2[4])+20):
                            if int(tupla[6])<=(int(tupla2[7])+20) or int(tupla[7])>=(int(tupla2[6])-20):
                                c=c+1
                                pass
                    if tupla2[12][0]=="-" and tupla2[12][3]=="-":
                        if int(tupla[4])<=(int(tupla2[5])+20) or int(tupla[5])>=(int(tupla2[4])-20):
                            if int(tupla[6])<=(int(tupla2[7])+20) or int(tupla[7])>=(int(tupla2[6])-20):
                                c=c+1
                                pass
                    if tupla2[12][0]=="-" and tupla2[12][3]=="+":
                        if int(tupla[4])<=(int(tupla2[5])+20) or int(tupla[5])>=(int(tupla2[4])-20):
                            if int(tupla[6])>=(int(tupla2[7])-20) or int(tupla[7])<=(int(tupla2[6])+20):
                                c=c+1
                                pass
                    if c==len(templist2):
                        templist2.append(tupla)
                        c=0                       
                        
                
                    #f3.write(str(tupla2)+"\n")
                    
                    
        for tupla2 in templist2:
            elements_for_identities=tupla2[8].split("/")
            identities=identities+int(elements_for_identities[0])
            elements_for_positives=tupla2[10].split("/")
            positives=positives+int(elements_for_positives[0])
            new_elem=elements_for_positives[1].split("(")
            alignment_len=alignment_len+int(new_elem[0])
            elements_for_gaps=tupla2[11].split("/")
            gaps=gaps+int(elements_for_gaps[0])
            score=score+float(tupla2[13])
            fragments=fragments+1
            query_length=int(tupla2[2])
            subject_length=int(tupla2[3])
            
            for element in tupla2:
                
                f2.write(element+"\t")
            f2.write("\n")
        if templist[0][0]!=ref:
            f3.write("# query"+"\t"+"subject"+"\t"+"query_length"+"\t"+"subject_length"+"\t"+"aln_legth"+"\t"+"identities"+"\t"+"positives"+"\t"+"gaps"+"\t"+"frame"+"\t"+"score"+"\t"+"n._of_fragments"+"\t"+"E-Value"+"\n")
        ref=templist[0][0]
        #print (query_length/3, alignment_len, subject_length/3)
        #if int(alignment_len)<=(query_length/3) and alignment_len<=(subject_length/3):    
        f3.write(templist[0][0]+"\t"+templist[0][1]+"\t"+str(round(query_length/3))+"\t"+str(round(subject_length/3))+"\t"+str(alignment_len)+"\t"+str(identities)+"\t"+str(positives)+"\t"+str(gaps)+"\t"+tupla2[12]+"\t"+str(round(score))+"\t"+str(fragments)+"\t"+reference_evalue+"\n")
        #elif int(alignment_len)>(query_length/3):
        #    f3.write(templist[0][0]+"\t"+templist[0][1]+"\t"+str(round(score))+"\t"+str(fragments)+"\t"+str(round(query_length/3))+"\t"+str(round(subject_length/3))+"\t"+str(round(query_length/3))+"\t"+str(identities)+"\t"+str(float(identities/query_length/3))+"\t"+str(positives)+"\t"+str(float(positives/query_length/3))+"\t"+str(gaps)+"\n")
        #elif int(alignment_len)>(subject_length/3):
        #    f3.write(templist[0][0]+"\t"+templist[0][1]+"\t"+str(round(score))+"\t"+str(fragments)+"\t"+str(round(query_length/3))+"\t"+str(round(subject_length/3))+"\t"+str(round(subject_length/3))+"\t"+str(identities)+"\t"+str(float(identities/subject_length/3))+"\t"+str(positives)+"\t"+str(float(positives/subject_length/3))+"\t"+str(gaps)+"\n")
        score=0
        fragments=0
        identities=0
        gaps=0
        positives=0
        alignment_len=0
        templist2=[]
        templist=[]
f1.close()
f2.close()
f3.close()

f1=open(".transcriptologsTEMP\\output2_part2.txt")
#f2=open("C:\Documents and Settings\Luca\Desktop\\provaBEST_HIT3.txt", "w")
f2=open(".transcriptologsTEMP\\output2_part3.txt", "w")

templist=[]
templist2=[]

for line in f1:
    if not line.startswith("#"):
        templist.append(line)
        
    elif line.startswith("#") and len(templist)!=0:
        for item in templist:
            elements=item.split()
            score=float(elements[9])
                
                #f2.write(elements[0]+"\t"+elements[1]+"\t"+elements[10]+"\t"+str(round(totalscore))+"\n")
            templist2.append((round(score),elements[0],elements[1],elements[2],elements[3],elements[4],elements[5],elements[6],elements[7],elements[8],elements[10],elements[11]))
                
        a=sorted(templist2,reverse=True)
                
        reference=a[0]
                
        ref_score=float(reference[0])
                
        for element in templist2:
                    
            current_score=element[0]
            if current_score >= ref_score:
                    
                f2.write(element[1]+"\t"+element[2]+"\t"+str(element[3])+"\t"+str(element[4])+"\t"+str(element[5])+"\t"+str(element[6])+"\t"+str(element[7])+"\t"+str(element[8])+"\t"+str(element[9])+"\t"+str(element[0])+"\t"+str(element[10])+"\t"+str(element[11])+"\n")
                    
        templist=[]
        templist2=[]
f1.close()
f2.close()

f1=open(".transcriptologsTEMP\\output2_part3.txt")

f2=open(".transcriptologsTEMP\\output2_part3_TEMP.txt", "w")

ALL=f1.readlines()
ALL.sort()
for line in ALL:
    f2.write(line)
f2.close()
f1.close()
f2=open(".transcriptologsTEMP\\output2_part3_TEMP.txt")
f3=open(".transcriptologsTEMP\\output2_part4.txt", "w")

p1="vuoto"
p2="vuoto"
p3=100000
previousline="primalinea"
templist=[]
for line in f2:
    elements=tuple(line.split())
    n1=elements[0]
    
    n2=elements[1]
    n3=elements[2]
    if len(templist)==0:
        templist.append(elements)
    
    elif not n1==p1:
        
        templist.sort(key=lambda obj:int(obj[2]), reverse=True)
        reference=templist[0][2]

        
        for item in templist:
            
            if int(item[2])>=int(reference):
                f3.write(item[0]+"\t"+item[1]+"\t"+item[2]+"\t"+item[3]+"\t"+item[4]+"\t"+item[5]+"\t"+item[6]+"\t"+item[7]+"\t"+item[8]+"\t"+item[9]+"\t"+item[10]+"\t"+item[11]+"\n")
        templist=[]
        templist.append(elements)
    
    elif n1==p1 and n2==p2 and n3!=p3:
        templist.append(elements)
    elif n1==p1 and n2!=p2:
        templist.append(elements)
        
        
    p1=elements[0]
    p2=elements[1]
    p3=elements[2]
f1.close()
f2.close()
f3.close()

f1=open(".transcriptologsTEMP\\output2_part4.txt")

f2=open(".transcriptologsTEMP\\output2_part5.txt", "w")

ALL=[]

for line in f1:
    elements=line.split()
    newelements=elements[1]+"\t"+elements[0]+"\t"+elements[2]+"\t"+elements[3]+"\t"+elements[4]+"\t"+elements[5]+"\t"+elements[6]+"\t"+elements[7]+"\t"+elements[8]+"\t"+elements[9]+"\t"+elements[10]+"\t"+elements[11]+"\n"
    ALL.append(newelements)

ALL.sort()
for element in ALL:
    f2.write(element)
f1.close()
f2.close()

f1=open(".transcriptologsTEMP\\output1_part5.txt")
f2=open(".transcriptologsTEMP\\output2_part5.txt")
f3=open(args.outputfile, "w")

sorghum=f1.readlines()
f1.close()
arab=f2.readlines()
f2.close()
newindex=0
f3.write("QUERY"+"\t"+"SUBJECT"+"\t"+"QUERY_LEN(1)"+"\t"+"SUBJECT_LEN(1)"+"\t"+"ALN_LEN(1)"+"\t"+"IDENTITIES(1)"+"\t"+"POSITIVES(1)"+"\t"+"GAPS(1)"+"\t"+"FRAME(1)"+"\t"+
         "QUERY_LEN(2)"+"\t"+"SUBJECT_LEN(2)"+"\t"+"ALN_LEN(2)"+"\t"+"IDENTITIES(2)"+"\t"+"POSITIVES(2)"+"\t"+"GAPS(2)"+"\t"+"FRAME(2)"+"\t"+"SCORE(1)"+"\t"
         "E-VALUE(1)"+"\t"+"FRAGMENTS(1)"+"\t"+"SCORE(2)"+"\t"+"E-VALUE(2)"+"\t"+"FRAGMENTS(2)"+"\n")
for relat in sorghum:
    fields=relat.split()
    s1=fields[0]
    s2=fields[1]
    s3=fields[2]
    s4=fields[3]
    s5=fields[4]
    s6=fields[5]
    s7=fields[6]
    s8=fields[7]
    s9=fields[8]
    s10=fields[9]
    s11=fields[10]
    s12=fields[11]
    for relat2 in arab[newindex:]:
        fields2=relat2.split()
        a1=fields2[0]
        a2=fields2[1]
        a3=fields2[2]
        a4=fields2[3]
        a5=fields2[4]
        a6=fields2[5]
        a7=fields2[6]
        a8=fields2[7]
        a9=fields2[8]
        a10=fields2[9]
        a11=fields2[10]
        a12=fields2[11]
        if s1==a1 and s2==a2:
            newindex=arab.index(relat2)
            f3.write(s1+"\t"+s2+"\t"+s3+"\t"+s4+"\t"+s5+"\t"+s6+"\t"+s7+"\t"+s8+"\t"+s9+"\t"+a3+"\t"+a4+"\t"+a5+"\t"+a6+"\t"+a7+"\t"+a8+"\t"+a9+"\t"+s10+"\t"+s12+"\t"+s11+"\t"+a10+"\t"+a12+"\t"+a11+"\n")
            break
        if s1<a1:
            break    
f1.close()
f2.close()
f3.close()

shutil.rmtree(".transcriptologsTEMP")

h=time.time()
print(str(h-g)+" seconds")