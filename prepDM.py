#! /usr/bin/python


DMname = 'weighted_unifrac_dm.csv'
mapfilename = 'newdeermap.csv'

mapfile = open(mapfilename, 'r')

linenum = 0

Dict = {}

for line in mapfile:
	if linenum > 0:
		line = line.strip('\n')
		ElementList = line.split('\t')
		Dict[ElementList[0]] = ElementList[16] + "_" + ElementList[47]
		print Dict
	linenum += 1
mapfile.close()


DM = open("weighted_unifrac_dm.csv", 'r')
line2 = 0
cellnum = 0
newlist = []	
for line in DM:
	if line2 == 0:
		line = line.strip('\n')
		DMList = line.split(',')
		for cell in DMList:
			if cellnum == len(DMList):
				newlist.append('\n')
				break
			else:
				if DMList[cellnum] in Dict:
					newlist.append(Dict[DMList[cellnum]])	
					cellnum += 1
				else:
					newlist.append(DMList[cellnum])
					cellnum += 1
		line2 += 1
	elif line2 > 0:
		break
DM.close()

String = ""
linenum = 0
DMagain = open("weighted_unifrac_dm.csv", 'r')
for line in DMagain:
	if linenum == 0:
		linenum += 1
	else:		
		line = line.strip('\n')
		LList = line.split(',')
		#del LList[0]
		#LList.insert(0, newlist[linenum])
		   # inserts string at beginning of list instead of separating each letter
		editedline = ','.join(LList)
		String = String + editedline + '\n'
		linenum += 1				
DMagain.close()

header = ','.join(newlist)
outfile = open("indiweight.csv", 'w')
outfile.write(header + '\n' + String)
outfile.close()

