####################################################################################
#
#	Propagate gVCF or VCF formati mutation results throught the gene network(Network-based Stratification)
#
#	Author:	Chengcheng Ma
#
#	Date:	2017-3-5
#	
#	Need Packages: 	Python 	getopt
#					R	org.HS.eg.db, data.table, getopt		
#####################################################################################



import os, getopt, sys

def usage():
	print '''
		This is a programm for doing network propagation from gVCF or VCF format data through gene networks
		such as: STRING, HumanNet, or other networks. 

		The format of the networks should be like(-m option):
		gene1		gene2		interaction_score
		or(-c option):
				gene1	gene2
		gene1	score	socre
		gene2	socre	score
	
		example: 	python NetworkPropagation.py -i VCF_path/gVCF_file -o output_dir -n network_form1
				python NetworkPropagation.py -i VCF_path/gVCF_file -o output_dir -m network_form2
	
		-i		input		"character"		input path of VCF files or gVCF file[must]
		-o		output		"character"		network propagated output dir[must]
		-n		network		"character"		network to be propagated[must -n or -m]
		-m		network2	"character"		the gene interaction matrix format[must -n or -m]
		-g		germline	"germline"		1 or 0 whether sparate the homozugous and heterzygous mutations[default not sepearate 0]
		-s		step_length	"numeroc		the step length of propagation[default 0.5]
		-c		contraction	"numeric"		the contraction limit[default 1e-6]	
	'''


opts, args = getopt.getopt(sys.argv[1:], "h:i:n:m:o:g:", ["help", "input", "network", "network2", "output",  "germline"])

for opt, value in opts:
	if opt in "-h":
		usage()
		sys.exit()
	if opt in  "-i":
		input = value
	if opt in  "-o":
		outdir = value
	if opt in "-n":
		network = value
	if opt in "-m":
		network2 = value
	if opt in "-g":
		germline = value
	if opt in "-s":
		step = value
	if opt in "-c":
		contract = value

if not("input" in locals().keys() and "outdir" in locals().keys()):
	usage()
	exit(1)

if not("network" in locals().keys() or "network2" in locals().keys()):
	usage()
	exit(1)
if "network" in locals().keys():
	if not(os.path.exists(input) and (os.path.exists(network))):
		print "Input file or network file not exist!"
		exit(1)
else:
	if not(os.path.exists(input) and (os.path.exists(network2))):
                print "Input file or network file not exist!"
                exit(1)

if os.path.isfile(input):
	gvcf = 1
	prefix = os.path.basename(input).split(".")[0]
elif os.path.isdir(input):
	gvcf = 0
	prefix = os.path.basename(os.path.normpath(input))
#print prefix

if not os.path.exists(outdir):
	os.mkdir(outdir)


if not "gvcf" in  locals().keys():
	gvcf = 0
if not "germline" in  locals().keys():
    germline = 0
if not "step" in  locals().keys():
    step = 0.5
if not "contract" in  locals().keys():
    contract = 1e-6

path = os.path.dirname(os.path.realpath(__file__))
#print path
if "network" in locals().keys():
	netbase = os.path.basename(network).split(".")[0]


if "network" in locals().keys():
	cmd1 = "Rscript %s/IntMatrix.R -i %s -o %s/%s.NetMatrix " % (path, network, outdir, netbase)
	#print cmd1
	os.system(cmd1)
	cmd2 = "Rscript %s/Vcf2Matrix.R -i %s -n %s/%s.NetMatrix -o %s/%s.Matrix -g %d -m %d" % (path, input, outdir, netbase, outdir, prefix, int(gvcf), int(germline))
	#print cmd2
	os.system(cmd2)
	cmd3 = "Rscript %s/NetProp.R -i %s/%s.Matrix -n %s/%s.NetMatrix -o %s/%s.Propagated" % (path, outdir, prefix, outdir, netbase, outdir, prefix)
	#print cmd3
	os.system(cmd3)
if "network2" in locals().keys():
	cmd1 = "Rscript %s/Vcf2Matrix.R -i %s -n %s -o %s/%s.Matrix -g %d -m %d" % (path, input,  network2, outdir, prefix, int(gvcf), int(germline))
	#print cmd1
	os.system(cmd1)
	cmd2 = "Rscript %s/NetProp.R -i %s/%s.Matrix -n %s -o %s/%s.Propagated" % (path, outdir, prefix,  network2, outdir, prefix)
	#print cmd2
	os.system(cmd2)


