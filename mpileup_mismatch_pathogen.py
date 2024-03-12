import sys
import random
import math
from math import sqrt


mindepth=float(sys.argv[1])
block_size=int(sys.argv[2])

nucleotides=['A','C','T','G']
def basefun(ref,inp):
	string =''
	for thebase in inp:
		if thebase.isalpha() and thebase.upper() in nucleotides:
			string += thebase.upper()
		else:
			if thebase in ['.',',']:
				string += ref
	return string

t_list=[]
n_list=[]
choice_list=[]

currentchrom=0
mmatch=0
mismatch=0
blockmatch=0
blockmismatch=0
for line in sys.stdin:
	newline =''
	col=line.split()
	chromosome=col[0]

	if len(col) <8: continue
	position=int(col[1])
	refbase=col[2]
	depth=col[3]
	depth2=col[6]
	if '*' in depth or '*' in depth2: continue
	#if refbase =='N' or int(depth) <1 or int(depth2) <1:continue
	if int(depth) < mindepth or int(depth2)<mindepth: continue
	pileup=basefun(refbase,col[4])
	pileup2=basefun(refbase,col[7])
	
	if len(pileup) < mindepth:continue
	if len(pileup2) < mindepth:continue
	if len(set(pileup))>1:continue
	if len(set(pileup2))>1:continue
	
	fullpileup=pileup+pileup2
	if 'C' in fullpileup and 'T' in fullpileup:continue
	if 'G' in fullpileup and 'A' in fullpileup:continue	

	if pileup[0] != pileup2[0]:
		t=1.0
		print chromosome,position,pileup,pileup2
		#print '1',line.rstrip('\n')
	elif pileup[0] == pileup2[0]:
		t=0.0
		#print '0',line.rstrip('\n')
	
	n=1.0
	#print pileup,pileup2,t
	t_list.append(t)
	n_list.append(n)	
	choice_list.append((chromosome,position,t,n))


########################################################

sumt_main=sum(t_list)
num_sites=sum(n_list)

mismatch_rate=sumt_main/num_sites

SE=math.sqrt((mismatch_rate*(1.0-mismatch_rate))/num_sites)

print 'Results:\t',mismatch_rate,'\t',SE,'\t',int(sumt_main),'\t',int(num_sites)

########################################################

lenmain=len(t_list)
sumt_main=sum(t_list)
sumn_main=sum(n_list)

num_SNPs=len(t_list)

try:
	Dp_main = sum(t_list) / sum(n_list)
except ZeroDivisionError:
	
	print 'Results:\t',mismatch_rate,'\t',SE,'\t',int(sumt_main),'\t',int(num_sites)
	exit(0)

if num_SNPs <50:
	print 'NA','\t','NA','\t',num_SNPs
	exit(0)
	

block_size=10000

mjlist=[]
t_list_b= []
n_list_b= []

t_list_copy=t_list
n_list_copy=n_list
jackknife_Dp=[]
jackknife_D=[]
jackknife_E=[]
current_chromosome = 0
blockcount = 0

bcounter=0
for line in choice_list:
	col=line
	chromosome = col[0]
	position = int(col[1])

	t = col[2]
	n = col[3]
	bcounter+=1
	if bcounter==1:
		prev_position=position
		#current_chromosome=chromosome
		blockstartindex=bcounter
	t_list_b.append(t)
	n_list_b.append(n)

	if False:
		continue
	else: #basepair blocks
		if chromosome != current_chromosome:
		
			if True:
				sumt=sumt_main - sum(t_list_b)
				sumn=sumn_main- sum(n_list_b)
				
				mjlist.append(len(t_list_b))
		                       
				D_pop =sumt / sumn
				jackknife_Dp.append(D_pop)
			t_list_b=[]
			n_list_b =[]
			current_chromosome = chromosome
			prev_position= 0
			continue
		elif position > (prev_position+block_size):
			

			
			sumt=sumt_main- sum(t_list_b)
			sumn=sumn_main- sum(n_list_b)

			D_pop =sumt / sumn
			mjlist.append(len(t_list_b))
		                       

			jackknife_Dp.append(D_pop)

			t_list_b=[]
			n_list_b =[]
			prev_position=position

noweighting=False

if noweighting:    
	n_groups=float(len(jackknife_Dp))      
	pseudovalues=[]
	for j in jackknife_Dp:
		pseudovalues.append( float((n_groups*Dp_main - (n_groups -1.0)*j)) )
	jackknife_estimate = float(sum(pseudovalues)) / n_groups

	sum_list=[]
	for f in pseudovalues:
		sum_list.append(float(f-jackknife_estimate)**2)

	Dp_SE= sqrt(  sum( sum_list ) / (n_groups*(n_groups-1.0))  )            
	
else:
	n_groups=float(len(jackknife_Dp))
	
	pseudovalues=[]
	for m,Dj in zip(mjlist,jackknife_Dp):
		pseudovalues.append(  ((num_SNPs-m) * Dj)/num_SNPs )
	jackknife_estimate = n_groups*Dp_main - float(sum(pseudovalues))
	
	#print jackknife_estimate
	jvallist=[]
	for mj,Dj in zip(mjlist,jackknife_Dp):
		hj=num_SNPs/mj
		thetaminusj=Dj
		tj=hj*Dp_main-(hj-1.0)*thetaminusj
		
		jval=((tj-jackknife_estimate)**2) / (hj-1.0)
		jvallist.append(jval)
	Dp_SE=sum(jvallist)/n_groups
	Dp_SE= sqrt(Dp_SE )


print 'Results:\t',mismatch_rate,'\t',Dp_SE,'\t',int(sumt_main),'\t',int(num_sites)


