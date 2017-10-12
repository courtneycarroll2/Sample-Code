#! /usr/bin/python
## Shows read counts from OTU table for different sample categories
import sys

def taxonomy(querystring, stringlist, begin2, end2, somesum, level, count):
	query = querystring.split('\t')
	tax_list = query[begin2:end2]
	tax_string = '\t'.join(tax_list)
	tax_string = tax_string.strip()
	#tax_assignment = query[begin2:end2]
	#tax_assignment = '\t'.join(tax_assignment)	
	#tax_assignment= tax_assignment.strip()
	#somekey = str(query[0]) + '\n' + tax_assignment
	if tax_string not in stringlist:
		stringlist.append(tax_string)
		level[tax_string] = (somesum)
		count += 1
		return level
	else:
		oldkey = level[tax_string]
		level[tax_string] = str(int(oldkey) + int (somesum))
		return level

def OTUtaxonomy(querystring, stringlist, begin2, end2, somesum, level, count):
	query = querystring.split('\t')	
	tax_list = query[:end2]
	tax_list2 = query[begin2:end2]
	tax_string = '\t'.join(tax_list2)
	tax_string = tax_string.strip()
	tax_list[1:begin2] = []
	tax_assignment = '\t'.join(tax_list)	
	tax_assignment= tax_assignment.strip()
	#somekey = str(query[0]) + '\n' + tax_assignment
	if tax_string not in stringlist:
		stringlist.append(tax_string)
		level[tax_assignment] = (somesum)
		count += 1
		return level
	else:
		oldkey = level[tax_assignment]
		level[tax_assignment] = str(int(oldkey) + int (somesum))
		return level



InFile = open(sys.argv[1], 'r')
OutFile = open(sys.argv[2], 'w')

linecount = 1
if sys.argv[3] == "AccuprimeMoBioFeces":
	begin1 = int(1)
	end1 = int(6)
elif sys.argv[3] == "AccuprimeMoBioC1":
	begin1 = 6
	end1 = 11
elif sys.argv[3] == "AccuprimeSaltingOutFeces":
	begin1 = 11
	end1 = 16
elif sys.argv[3] == "AccuprimeSaltingOutC1":
	begin1 = 16
	end1 = 21
elif sys.argv[3] == "AccuprimeZymoFeces":
	begin1 = 21
	end1 = 25
elif sys.argv[3] == "AccuprimeZymoC1":
	begin1 = 25
	end1 = 30
elif sys.argv[3] == "FivePrimeMoBioC1":
	begin1 = 30
	end1 = 35
elif sys.argv[3] == "FivePrimeSaltingOutC1":
	begin1 = 35
	end1 = 40
elif sys.argv[3] == "FivePrimeZymoFeces":
	begin1 = 40
	end1 = 45
elif sys.argv[3] == "FivePrimeZymoC1":
	begin1 = 45
	end1 = 50

OTUs = {}
OTUcount = 0
kingdom = {}
kingdomcount = 0
phylum = {}
phylumcount = 0
Class = {}
classcount = 0
order = {}
ordercount = 0
family = {}
familycount = 0
genus = {}
genuscount = 0
species = {}
speciescount = 0
kstringlist = []
pstringlist = []
cstringlist = []
ostringlist = []
fstringlist = []
gstringlist = []
sstringlist = []
lstringlist = []



for line in InFile:
	if linecount == 18336:
		break

	elif line[0] != "#":
		line = line.strip()
		MyList = line.split(',')
		for cell in MyList:
			if cell == "":
				cell = "NA"
		MyString = '\t'.join(MyList)
		MyString = MyString.strip()
		SpecList = MyList[begin1:end1]
		SpecList2 = [int(x) for x in SpecList]
		reads = sum(SpecList2) 
		#reads = int(MyList[1]) + int(MyList[2]) + int(MyList[3]) + int (MyList[4]) + int(MyList[5])
		#reads = reduce(lambda q,p: p+q, SpecList)
		if reads != 0:
			k = taxonomy(MyString, kstringlist, 50, 51, reads, kingdom, kingdomcount)
			p = taxonomy(MyString, pstringlist, 50, 52, reads, phylum, phylumcount)
			c = taxonomy(MyString, cstringlist, 50, 53, reads, Class, classcount)
			o = taxonomy(MyString, ostringlist, 50, 54, reads, order, ordercount)
			f = taxonomy(MyString, fstringlist, 50, 55, reads, family, familycount)
			g = taxonomy(MyString, gstringlist, 50, 56, reads, genus, genuscount)
			s = taxonomy(MyString, sstringlist, 50, 57, reads, species, speciescount)	
			l = OTUtaxonomy(MyString, lstringlist, 50, 58, reads, OTUs, OTUcount)

	linecount += 1
 			
OutFile.write(sys.argv[3] + "stats:\nKingdom\n")
for key in kingdom:
	OutFile.write(key + "\t\t\t" + str(kingdom[key]) + '\n')

OutFile.write('\n\nPhylum\n')
for key in phylum:
	OutFile.write(key + "\t\t\t" + str(phylum[key]) + '\n')

OutFile.write('\n\nClass\n')
for key in Class:
	OutFile.write(key + "\t\t\t" + str(Class[key]) + '\n')

OutFile.write('\n\nOrder\n')
for key in order:
	OutFile.write(key + "\t\t\t" + str(order[key]) + '\n')

OutFile.write('\n\nFamily\n')
for key in family:
	OutFile.write(key + "\t\t\t" + str(family[key]) + '\n')

OutFile.write('\n\nGenus\n')
for key in genus:
	OutFile.write(key + "\t\t\t" + str(genus[key]) + '\n')

OutFile.write('\n\nSpecies\n')
for key in species:
	OutFile.write(key + "\t\t\t" + str(species[key]) + '\n')

OutFile.write('\n\nOTUs\n')
for key in OTUs:
	OutFile.write(key + "\t\t\t" + str(OTUs[key]) + '\n')

InFile.close()
OutFile.close()

