import math
from Bio import SeqIO

def euclid(dic,dict1):
	x=[]
	y=[]
	for i in dic.values():
		x.append(i)

	for j in dict1.values():
		y.append(j)

	n=len(x)
	dist=0
	for i in range(n):
		dist=dist+((x[i]-y[i])**2)

	dist=math.sqrt(dist)

	return dist

def word(a,b):
	l=[]
	k=3
	n=len(a)-k
	i=0
	while i<=n:
		l.append(a[i:i+k])
		i=i+1


	w=[]
	m=len(b)-k
	j=0
	while j<=m:
		w.append(b[j:j+k])
		j=j+1


	e=set(l)
	r=set(w)
	h=e.union(r)

	dic={}
	for each in h:
		dic[each]=0
		for i in l:
			if i==each:
				if i not in dic.keys():
					dic[i]=1
				else:
					dic[i]=dic[i]+1

	
	dict1={}
	for each in h:
		dict1[each]=0
		for i in w:
			if i==each:
				if i not in dic.keys():
					dict1[i]=1
				else:
					dict1[i]=dict1[i]+1


	val=euclid(dic,dict1)

	return val<30

def readfile():

	seq=[]
	ide=[]
	for seq_record in SeqIO.parse("sup.fasta","fasta"):
		seq.append(str(seq_record.seq))
		ide.append(seq_record.id)

	print(len(seq))
	return seq,ide


def main():
	seq,ide=readfile()

	rep={}
	vis=[]
	n=len(seq)

	vis=[False for i in range(n)]

	for i in range(n):
		if vis[i]==False:
			for j in range(i+1,n):
				if vis[j]==False and word(seq[i],seq[j]):
					vis[i]=True
					vis[j]=True
					if i not in rep.keys():
						rep[i]=[j]
					else:
						rep[i].append(j)

	for i in range(n):
		if vis[i]==False:
			rep[i]=i

	finalrep=set()
	for i in rep.keys():
		finalrep.add(ide[i])

	print(finalrep)



if __name__ == '__main__':
	main()

