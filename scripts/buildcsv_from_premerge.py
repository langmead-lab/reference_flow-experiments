import glob, os, sys, subprocess

dataDir=sys.argv[1]
outPrefix=sys.argv[2]
#thresh=[5]
thresh=[5,10,15,20,25,30]
#outName=sys.argv[2]
#suffix=".fa.msh"

total_recall=0.0

samples=glob.glob("{}/*t{}-*.log".format(dataDir,thresh[0]))

fName=[ i.split("_pass2-t{}-".format(thresh[0]))[0] for i in samples ]

#print (glob.glob("{}/*t{}*_accuracy.txt".format(dataDir,thresh)))
for i,t in enumerate(thresh):
	outFile=outPrefix+'_t{}.csv'.format(t)
	outList=[]
	for j,f in enumerate(fName):

		outName = f.split("/")[-1]		

		p1_time=0.0
		p2_time=0.0
		p1_log = "{}.log".format(f)
		p2_log = "{}_pass2-t{}-1s0_b1.log".format(f,t)

		with open(p1_log) as myFile:
			for line in myFile:
				if line:
					line=line.strip()
					if "Elapsed" in line:
						time_string=line.split(" ")[-1]
						time_array=time_string.split(":")
						if len(time_array)==2:
							p1_time = float(time_array[0])*60.0+float(time_array[1])
						elif len(time_array)==3:
							p2_time = float(time_array[0])*60*60+float(time_array[1])*60+float(time_array[2])

		with open(p2_log) as myFile:
			for line in myFile:
				if line:
					line=line.strip()
					if "Elapsed" in line:
						time_string=line.split(" ")[-1]
						time_array=time_string.split(":")
						if len(time_array)==2:
							p2_time += float(time_array[0])*60.0+float(time_array[1])
						elif len(time_array)==3:
							p2_time += float(time_array[0])*60*60+float(time_array[1])*60+float(time_array[2])

		#print("{} {}".format(p1_time,p2_time))

		p1_accuracy = "{}_premerge_acc.txt".format(f)
		p1_tp = 0.0
		p1_sall = 0.0
		with open(p1_accuracy) as myFile:
			for line in myFile:
				if line:
					line=line.strip()
					if 'sensitivity_all' in line:
						line = line.split()
						p1_tp = int(line[3][1:])
						p1_sall = int(line[5][:-1])

		p2_regex = "{}_pass2-t{}-*premerge_acc.txt".format(f,t)
		p2_tp = 0.0
		p2_sall = 0.0
		for k in (glob.glob(p2_regex)):
			with open(k) as myFile:
				for line in myFile:
					if line:
						line=line.strip()
						if 'sensitivity_all' in line:
							line = line.split()
							p2_tp+=int(line[3][1:])
							p2_sall+=int(line[5][:-1])

		outList.append("{}, {}, {}, {}, {}, {}, {}\n".format(outName, p1_time, p2_time,p1_tp, p1_sall, p2_tp,p2_sall))
			
	with open(outFile,'w') as myFile:
		for l in outList:
			myFile.write(l)
