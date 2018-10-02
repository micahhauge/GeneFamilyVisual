#!/usr/bin/python3

# Authors: Anthony Davis, Hyrum D. Carroll
# Description: Parses FASTA files (in UCSC Genome Browser CDS format) into JSON records
# Note: execute as python3 -o ....py to avoid debugging code
# To-Do:
#        Make start and stop reading frames {^, <, >} be globals
#        Do we need both start and stop frames?
#        Determine if were running from the terminal or a web form

import os, string, sys, cgi
import binascii
#import binascii, mkdir
from Bio import SeqIO  # requires BioPython
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline   # requires Clustal Omega (with the executable as clustalo)
import argparse #import argaparse library for addArgv, --afa
import pprint  # for debugging only

UTR_EXON = -2
MIN_AA_LEN = 100   # minimum amino acid length to consider

UPLOADED_FILES_DIR = 'uploadedFiles/'

# For DEBUGGING in cgi
import cgitb;


#-------------------------------------------------

CLUSTAL_INPUT_BASE_FILENAME = "clustalInput.fa"

def debug(msg):
    if __debug__:
        print(f"<br/>\nDEBUG: {msg}<br/>\n")

#creating an object to hold all the information for each gene
class Gene(object):
    def __init__(self, name):
        #name of gene
        self.name = name

        # NOT USED? 20180906
        # # array that will hold exon lengths 
        # self.ExonLens = []

        # length of AA (unaligned)
        self.AA_len = 0

        # string to store the translated (unaligned) gene
        self.AA = ''

        #string to store the CDS region and lengths
        self.CDS = ''
        self.CDS_len = []

        # string to store the translated (aligned) gene
        self.Aligned_str = ''

        #string to store the frame string
        self.Frame_str = ''

        #array to keep track of exon frames
        self.Frames = []

        #line color for this gene
        self.color = ''
        
# dictionary of codon to amino acid
dnaToProtD = {
    "TTT":'F', "TTC":'F', "TTA":'L', "TTG":'L',
    "CTT":'L', "CTC":'L', "CTA":'L', "CTG":'L',
    "ATT":'I', "ATC":'I', "ATA":'I', "ATG":'M',
    "GTT":'V', "GTC":'V', "GTA":'V', "GTG":'V',

    "TCT":'S', "TCC":'S', "TCA":'S', "TCG":'S',
    "CCT":'P', "CCC":'P', "CCA":'P', "CCG":'P',
    "ACT":'T', "ACC":'T', "ACA":'T', "ACG":'T',
    "GCT":'A', "GCC":'A', "GCA":'A', "GCG":'A',

    "TAT":'Y', "TAC":'Y', "TAA":'', "TAG":'',
    "CAT":'H', "CAC":'H', "CAA":'Q', "CAG":'Q',
    "AAT":'N', "AAC":'N', "AAA":'K', "AAG":'K',
    "GAT":'D', "GAC":'D', "GAA":'E', "GAG":'E',

    "TGT":'C', "TGC":'C', "TGA":'', "TGG":'W',
    "CGT":'R', "CGC":'R', "CGA":'R', "CGG":'R',
    "AGT":'S', "AGC":'S', "AGA":'R', "AGG":'R',
    "GGT":'G', "GGC":'G', "GGA":'G', "GGG":'G'
}

# display an error or warning message to HTML or the terminal
def alert(msg):
    # FUTURE: Determine if were running from the terminal or a web form
    print(msg)
    

# convert DNA into protein
# initialize aminos as a list of amino acids from each slice 3 letter DNA strand
# initialize protein as a string from join aminos list
# return protein
# Note: Assumes that all chracters are nucleotides (AGCT) and not gaps

def dna_to_prot(strand):
    if __debug__: print(f"dna_to_prot( (len: {len(strand)}) {strand})")
    if len( strand ) % 3 != 0:
        extraNucleotides = len( strand ) % 3
        print(f'ALERT: The length of the DNA string is NOT a multiple of 3, so the last {extraNucelotides} nucleotides ({strand[ -extraNucleotides:]}) will be ignored (DNA: {strand})')
        strand = strand[ :-extraNucleotides]
        
    aminos = [ dnaToProtD[ strand[i:i+3] ] for i in range(0, len(strand), 3) ]
    protein = "".join(aminos)
    return protein

# make a new directory with a random name (in the files directory)
# if files directory has not made, create one and set permission to 0o777
# initialize newDir as a random string, such as 424e069cb11ecc4f5712c3863c9cfba1
# initialize a path that is returned and used in Main()
# store a path of newDir into path string
# NOTE: This 0o777 is used for DEBUGGING purpose only, recommend use 0o755
# Because the server user permission is wwww-data, if using 0o755 is hard to remove test files
def mkDir():
    baseDir = UPLOADED_FILES_DIR
    if not os.path.exists(baseDir):
        os.makedirs(baseDir)
        os.chmod(path, 0o755)
        if __debug__: os.chmod(baseDir, 0o777) # DEBUGGING
        
    newDir = binascii.hexlify(os.urandom(16)).decode()
    path = baseDir + '/' + newDir
    if not os.path.exists(path):
        os.makedirs(path)
        os.chmod(path, 0o755)
        if __debug__: os.chmod(path, 0o777) # DEBUGGING

    return path

# This function will create a json file and potentially set permissions
def writeJson(dictJson, jsonFilename):
    try:
        jsonFileObj = open(jsonFilename, 'w')
        jsonFileObj.write(json.dumps(dictJson))
        jsonFileObj.close()
    except:
        alert("ERROR: Can not open and write to JSON file," + jsonFilename + "!")
        sys.exit(1)
        
    os.chmod(jsonFilename, 0o755)
    if __debug__: os.chmod(jsonFilename, 0o777)



def exonLocs2JSON( gene_dic ):
    debug(f"exonLocs2JSON( {gene_dic} )")
    if __debug__: pprint.pprint( gene_dic )
    
    # Make Final 2D array
    CombinedLists = []  # 1st index: genes; 2nd index items are: name, CDSExonCount, startFrame, stopFrame, CDSlength, exonLength, then AA alignment indicies
    geneIndex = 0
    # each gene is a row in CombinedLists
    # append each element in gene.Frame into CombinedLists
    for name, gene in gene_dic.items():
        CombinedLists.append([name])
        for col in gene.Frames:
            CombinedLists[geneIndex].append(col)
        geneIndex += 1

    # initialize preFinal to ...
    # for each list in CombinedLists as frameInfoList
    #   intialize name1 as first element in list which is name of gene
    #   append name into preFinal
    #   for each element in frameInfoList except first element name
    #       if alignmentIndex_firstExon is UTR_EXON, then append alignmentIndex_firstExon and exonSize
    #       else if count == 0
    #           for  in CombinedLists
        
    preFinal = []
    count = 0
    for frameInfoList in CombinedLists:
        name1 = frameInfoList[0]
        preFinal.append([name1])
        for col in frameInfoList[1:]:
            rowCount = count + 1
            rowOn = 0
            exonNum = col[0]  # coding exon
            startFrame = col[1]
            endFrame = col[2]
            # = col[3]
            exonSize = col[4] # coding length
            alignmentIndex_firstExon = col[5]
            
            if alignmentIndex_firstExon is UTR_EXON:  # if the first exon is only a UTR
                preFinal[count].append([alignmentIndex_firstExon, exonSize])
            elif count == 0:
                preFinal[count].append([alignmentIndex_firstExon, exonSize])
                for Nrow in CombinedLists[rowCount:]:
                    name2 = Nrow[0]
                    for Ncol in Nrow[1:]:
                        if startFrame == Ncol[1] and endFrame == Ncol[2]:  # start and end frames match
                            if (abs(exonSize - Ncol[4]) % 3) is 0: # lengths are the same ORF
                                if abs(alignmentIndex_firstExon - Ncol[5]) <= 10:
                                    # each of first alignment positions for each gene are within 10 of each other
                                    preFinal[count][exonNum+1].append([name2, Ncol[5], Ncol[4], Ncol[0]])
            else:
                found = False
                for row1 in preFinal:
                    if found is False:
                        lastName = row1[0]
                        colOn = 0
                        for col1 in row1[1:]:
                            if len(col1) > 2:
                                count3 = 0
                                while count3 < len(col1)-2 and found is False:
                                    if name1 == col1[count3+2][0] and exonNum == col1[count3+2][3]:
                                        found = True
                                        preFinal[count].append([-1, [rowOn, colOn, col1[0]]])
                                    count3 += 1
                            colOn += 1
                    rowOn += 1

                if found is False:
                    preFinal[count].append([alignmentIndex_firstExon, exonSize])
                    for Nrow in CombinedLists[rowCount:]:
                        name2 = Nrow[0]
                        for Ncol in Nrow[1:]:
                            if startFrame == Ncol[1] and endFrame == Ncol[2]:
                                if (abs( exonSize - Ncol[4]) % 3) is 0:
                                    if abs(alignmentIndex_firstExon - Ncol[5]) <= 10:
                                        preFinal[count][exonNum+1].append([name2, Ncol[5], Ncol[4], Ncol[0]])
            rowCount += 1
        count +=1



    # ###DEBUGGING SECTION TO SEE THE 2D ARRAY OF EXON LENGTHS
    # i = 1
    # x = 100
    # y = 20
    # temp = '<svg width=300% height=300%>\n'
    # 
    # #find out how many nonCDS exons there to find how far to shift SVG
    # maxStart = 0
    # for row in preFinal:
    #     startNum = 0
    #     for col in row[1:]:
    #         if col[0] == UTR_EXON:
    #             startNum += 1
    #     if startNum > maxStart:
    #         maxStart = startNum
    # # print('max start is', maxStart)
    # 
    # with open("createSVGtemp.html","w") as printSVG:
    #     printSVG.write("<html>\n")
    #     printSVG.write("<body>\n")
    #     printSVG.write("<p>dataParsing is completed, but still upgrade graph design</p>\n</body>\n</html>");


# # This function will make a html file to display result of code
# # It is not completed yet, depend on the javascript
# # NOTE: This 0o777 is used for DEBUGGING purpose only, recommend use 0o755
# # Because the server user permission is wwww-data, if using 0o755 is hard to remove test files
# def mkPage(path):
#     printPage = ("""<!DOCTYPE html>
# <html>
# <head>
#     <title>Result Upload File</title>
#     <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.0/jquery.min.js"></script>
#     <script type="text/javascript" src="../openJson.js"></script>
# </head>
# <body>
#     <p id="myJson">Thank you for using this upload. Unforunately, the back end is not finish yet...</p>
# </body>
# </html>"""%("jsonFile.json"))
#     newPage = path + "/result.html" 
#     open(newPage, 'w').write(printPage)
#     #os.chmod(newPage, 0o755)
#     os.chmod(newPage, 0o777)# DEBUGGING

# Processes all of the FASTA files in path (a directory)
# afa is an optional alignment filename or is False
def dataProcess( path, afa ):
    debug(f"dataProcess({path},{afa})")
    # create a dictionary to hold the Gene objects
    gene_dic = {}

    # loop through every file in the directory
    # For each file, initalize a gene object as follow;
    #    gene.name:      the name of file
    # #   gene.ExonLens: using SeqIo.parse to get exon and its length
    #    gene.CDS:       remove the UTR, by strip lowercase letter in exon
    #    gene.CDS_len:   append of each length of CDS after strip lowercase
    #    gene.Frames:    append of those follow:
    #        CDSExonCount: count of each ExonCount
    #        startFrame: if CDSExonCount is zero startFrame is '^'
    #                    else get the ORF from the previous exon
    #        stopFram:   if CDSLength % 3 == 0 --> '^'
    #                    if CDSLength % 3 == 1 --> '<'
    #                    if CDSLength % 3 == 2 --> '>'
    #        CDSlength:  if CDSExonCount == 0 --> len(CDS)
    #                    if CDSExonCount not 0 check prev-Frame is '<' --> -2
    #                    if CDSExonCount not 0 check prev-Frame is '>' --> -1
    #        exonLength: length of exon
    #    gene.CDS:       CDSExon with trim off 1-2 nucleotides if the coding region is not a multiple of 3
    #    gene.AA:        string of converted protein from called dna_to_prot(CDSExon)
    #    gene.AA_len:    length of AA string
    #  store each gene into gene_dic with its name
    
    # get the files in the directory
    listing = os.listdir(path)
    debug(f"listing: {listing}")
    
    clustalInputFilename = path + "/" + CLUSTAL_INPUT_BASE_FILENAME
    
    for infile in listing:
        debug(f"infile: {infile}")
        
        # if infile[-3] != ".fa":
        #     # not a fasta file extension (for command-line use)
        #     continue
        
        # open the file
        with open(os.path.join(path, infile), 'r') as genes:
            # initalize the counter to zero for each gene/file
            CDSExon = ''
            CDSExonCount = 0

            # get the name of the gene (from the filename) and create an empty object
            name = os.path.splitext(infile)[0]
            gene = Gene(name)


            # go through the file, get the gene name and then
            # get the length of all exons and put into array
            for seq_record in SeqIO.parse(genes, "fasta"):
                startFrame = ''
                stopFrame = ''

                # # add the exon length into ExonLens array
                # gene.ExonLens.append(len(seq_record))

                # convert Bio.SeqRecord to fasta format to get the sequence
                # in string format to check letter casing
                exon = seq_record.format("fasta")
                # this removes the first line of the fasta file and all newlines
                exon = exon.split('\n', 1)[-1].replace('\n', '')

                exonLength = len(exon)
                # remove the UTR (lowercase letters)
                CDS = exon.strip(string.ascii_lowercase)
                CDSExon += CDS
                CDSlength = len(CDS)
                gene.CDS_len.append(CDSlength)

                if CDSExonCount == 0 :
                    startFrame = '^'
                else:
                    startFrame = gene.Frames[CDSExonCount-1][2]  # get the ORF from the previous exon
                    # adjust the coding length to account for the reading frame from the previous exon
                    if startFrame == '<':
                        CDSlength -= 2
                    elif(startFrame == '>'):
                        CDSlength -= 1
                        
                # determine the reading frame
                if(CDSlength % 3 == 0):
                    stopFrame = '^'
                elif(CDSlength % 3 == 1):
                    stopFrame = '<'
                else:
                    stopFrame = '>'

                debug(f"CDSExonCount: {CDSExonCount}, startFrame: {startFrame}, stopFrame: {stopFrame}, CDSlength: {CDSlength}, exonLength: {exonLength}")
                gene.Frames.append([CDSExonCount, startFrame, stopFrame, CDSlength, exonLength])
                CDSExonCount += 1

            # trim off 1-2 nucleotides if the coding region is not a multiple of 3
            if len(CDSExon) % 3 != 0:
                excessNucleotides = len(CDSExon) % 3
                debug(f"Trimming {excessNucleotides} nucleotides from {CDSExon}")
                CDSExon = CDSExon[:-excessNucleotides]
                
            gene.CDS = CDSExon
            gene.AA = dna_to_prot(CDSExon)
            gene.AA_len = len(gene.AA)
            debug(f"len: {gene.AA_len}, AA: {gene.AA}")

            # Check to make sure all sequences are at least 100 amino acids long
            # Skip genes that are too short
            if gene.AA_len < MIN_AA_LEN:
                alert("Gene", name, " is too short with just", gene.AA_len, "amino acids (the minimum is" + str(MIN_AA_LEN) + ")!")
                continue
            
            if not afa:
                 # NOT USED? 20180906
                 # # create a directory to store AA fasta files
                 # aa_path = path + '/AA/'
                 # if not os.path.exists(aa_path):
                 #     os.makedirs(aa_path)
                 #     os.chmod(aa_path, 0o777) # DEBUGGING
                 
                 # FUTURE: open the file once, append to it in this loop, then close it once
                 # create one giant fasta file with all AA sequences to pass to clustal omega
                 with open( clustalInputFilename, "a") as ALL_AA:
                     ALL_AA.write('>' + name + '\n')
                     ALL_AA.write(gene.AA + '\n')
                     ALL_AA.close()

        gene_dic[name] = gene

    # if there is --afa on the command line, use the MSA from user
    # else use the alignment from the path in file folder
    # Either option: store Aligned_str into each gene in the gene_dictionary
    out_file = ""
    if not afa:
        out_file= path + "/aligned.fasta"
        # call clustal omega
        clustalo = ClustalOmegaCommandline( infile=clustalInputFilename, outfile=out_file, auto=True, cmd="clustalo" )
        clustalo()
    else:
        out_file = "./" + afa
    
    # read in the aligned sequences and store them in their respective objects
    with open(out_file, 'r') as clustal:
        for line in clustal:
            line = line.rstrip()
            if line[0] is '>':
                name = line[1:] # remove the > character at the beginning of the line
            else:
                if name not in gene_dic:
                    print(f'ERROR: Found {name} in {out_file} but that name was not in the input FASTA files!')
                    pprint.pprint(gene_dic)
                    sys.exit(1)
                gene_dic[name].Aligned_str += line
                
    # trim off 1-2 nucleotides if the coding region is not a multiple of 3
    gaplessLen = len( gene_dic[name].Aligned_str.replace("-","") )
    if gaplessLen % 3 != 0:
        excessNucleotides = gaplessLen % 3
        debug(f"Trimming {excessNucleotides} nucleotides from {gene_dic[name].Aligned_str}")
        gene_dic[name].Aligned_str = gene_dic[name].Aligned_str[:-excessNucleotides]


    if afa:
        # verify that the user-supplied alignment file matches exactly with the translated versions of the user-supplied input files
        # for each gene
        #     compare AA with Aligned_str (except for "-"s) (by making a temporary copy of Aligned_str that has all "-"s replaced with "")
        for name in gene_dic:
            if __debug__: print(f"Checking alignment of {name}...<br/>")
            gene = gene_dic[name]
            gaplessAligned = gene.Aligned_str.replace("-","")
            print( "Amino acids (len:", str(len(gene.AA)) + "):", gene.AA,"<br/>")
            print( "DNA (len:", str(len(gaplessAligned)) + "):", gaplessAligned, "(any gaps were removed for this error message)","<br/>")
            translated = dna_to_prot(gaplessAligned)
            
            if gene.AA != translated:
                print( "Translated version of input files for", name, "does not match with the supplied alignment:","<br/>")
                print( "Amino acids (len:", str(len(gene.AA)) + "):", gene.AA,"<br/>")
                print( "Alignment (translated) (len:", str(len(translated)) + "):", translated, "<br/>")
                print( "DNA (len:", str(len(gaplessAligned)) + "):", gaplessAligned, "(any gaps were removed for this error message)","<br/>")

                firstMismatch = -1
                for i in range( len(gene.AA) ):
                    if gene.AA[i] != translated[i]:
                        firstMismatch = i
                        break
                print( "First mismatch at index", firstMismatch)    
                sys.exit(1)    

                    
    # for each gene in the gene_dictionary
    #   intialize frame_len to get integer of each gene.CDS_len divide 3
        #   while frame_len is zero
        #       go to next exon then append UTR_EXON into gene.Frames
        #       also break if exon is more than length of gene.CDS_len
        #       else frame_len is next exon of gene.CDS_len divide 3
        #
        #   for each char in gene.Aligned_str
        #       if start is False, then set it to be true then append count into gene.Frame
        #       else if frame_len is zero, go to next exon then append count into gene.Frames and get frame_len from gene.CSD_len divide 3
        #       while frame_len is zero
        #           append UTR_EXON into gene.Frame then go to next exon
        #           get frame_len from gene.CSD_len divide 3
        #       else set frame_len to be -1
        #       increment count
        # 
        #   while exon is less than length of gene.CSD_len minus 1
        #       go to next exon then append UTR_EXON into gene.Frames
        
    for name, gene in gene_dic.items():
        debug(f"name: {name}, gene: {gene}")
        count = 0  # alignment index
        exon = 0   # index of the current        exon for this gene
        firstTime = True
        
        if len( gene.CDS_len) == 0:
            # no valid exons, go to the next gene
            continue
            
        # frame_len = int(gene.CDS_len[exon]/3)
        frame_len = gene.CDS_len[exon] // 3
        debug(f"frame_len: {frame_len}")

        # process all of the 5' exons that are entirely UTR (i.e., have at most 2 nucleotides)
        while frame_len == 0:  # only UTR
            gene.Frames[exon].append(UTR_EXON)
            debug(f"5' UTR exon (index {exon})")
            exon += 1
            if exon < len( gene.CDS_len):
                frame_len = gene.CDS_len[exon] // 3
            else:
                break

        # record the alignment index for
        # 1) the first AA and
        # 2) the last AA for each exon
        debug(f"Aligned_str: {gene.Aligned_str}")
        for char in gene.Aligned_str:
            if char is '-':
                pass
            else:
                if firstTime:
                    firstTime = False
                    # FUTURE: can this be pulled out of the for loop (i.e., can the frame_len == 0 for the first AA)?
                    gene.Frames[exon].append( count )  # record the index of the alignment
                    debug(f"exon index {exon}'s alignment starts at index {count}")
                elif frame_len is 0:  # done processing this exon
                    exon += 1
                    gene.Frames[exon].append( count )  # record the index of the alignment
                    debug(f"exon index {exon}'s alignment ends at index {count}")
                    frame_len = gene.CDS_len[exon] // 3
                    while frame_len == 0: # only UTR 
                        debug(f"3' UTR exon (index {exon})")
                        gene.Frames[exon].append(UTR_EXON)
                        exon += 1
                        frame_len = gene.CDS_len[exon] // 3
                else:
                    frame_len -= 1
            count += 1
           
        # process all of the 3' exons that are entirely UTR
        while exon < len(gene.CDS_len)-1:
            exon += 1
            gene.Frames[exon].append(UTR_EXON)
            debug(f"3'b UTR exon (index {exon})")

    exonLocs2JSON( gene_dic )

# Create a unique directory for this instance
# Get required FASTA input files from either CGI form or from the command-line
# Copy FASTA input files into unique directory (using the same file names)
# Process those input files
def main( form ):
    # os.system("chmod -R 777 " + UPLOADED_FILES_DIR)
    # os.system("echo 'DEBUGGING:' > /tmp/gfv.tmp; chmod 777 /tmp/gvf.tmp")

    # make a new unique directory in UPLOADED_FILES_DIR
    uniqueDir = mkDir()

    commandLine = False

    # if called from a CGI form
    # initialize fileItems as list of filename through HTML Form
    # for each fileItem in fileItems list do:
    # copy all fileItem into uniqueDir
        
    if 'GATEWAY_INTERFACE' in os.environ:  # Called from a CGI form
        fileItems = form['filename[]']
        for fileItem in fileItems:
            if fileItem.file:
                fn = os.path.basename(fileItem.filename.replace("\\", "/"))  # change Windows filenames
                try:
                    # copy the file
                    open(uniqueDir + '/' + fn, 'wb').write(fileItem.file.read())
                except:
                    # display error message
                    print("""Content-Type: text/html\n\n
                          <html>
                          <body>
                              <p>%s</p>
                          </body>
                          </html>
                          """ %(uniqueDir+'/'+fn))
                    sys.exit()
                                        
        
    else:         # called from a command-line
            # set commandLine to be True
            # initialize list fileItems from args.files which user provide
            # for each fileItem in fileItems list do:
            # copy all fileItem into uniqueDir
        
        # get file(s) from the command-line
            #
        # get input file(s) from either the CGI form or from the command-line
        #

        # initialize args to hold all arugment user provide
        parser = argparse.ArgumentParser()
        # set aligned file as an optional argument for user
        parser.add_argument("-a", "--afa", metavar="", help="Use alignment file instead of Clustal Omega alignment")
        # set list of fasta file to run program
        parser.add_argument("files",nargs="*", help="fasta file(s)")
        args = parser.parse_args()

        #pprint.pprint( args )
        
        # program require a list of fasta file to run a program
        # so check if user provide fasta file 
        # if not, print usage for user and terminate program

        if not args.files:
            print("usage: dataParsing.py [-h]/[--help]")
            print("usage: dataParsing.py [UCSC Genome Browser CDS fasta files]")
            print("usage: dataParsing.py [-a]/[--afa] [aligned file] [UCSC Genome Browswer CDS fasta files]")
            sys.exit(0)
            
        commandLine = True
        fileItems = args.files
        for fileItem in fileItems:
            fn = os.path.basename(fileItem)
            try:
                op = open(fileItem,'r').read()
                
                # Need to change 'wb' to 'w' due to this error message
                # 'str' does not support the buffer interface
                open(uniqueDir + '/' + fn, 'w').write(op)
            except:
                sys.exit()
    
        # initialize afa as filename of aligned file
        # if user do NOT provide a file, afa will set as False value
             
        afa = args.afa      #get the aligned provide, or set it as False if user dont provide the aligned
        
    dataProcess( uniqueDir, afa )

    # initialize redirectStr to display output of createSVGtemp.html if called from a CGI form
        
    redirectStr="""Content-Type: text/html\n\n
<html>
    <head>
        <meta http-equiv="refresh" content="0; url=createSVGtemp.html" />
    </head>
    <body>
        <p>Redirecting to <a href="createSVGtemp.html">createSVGtemp.html</a></p>
    </body>
</html>
"""
    
    # if called from CGI form use redirectStr
    # else print meassge of complete program
    if not commandLine:
        print(redirectStr)
    else:
        print("Successfully completed")

if __name__ == "__main__":

    form = cgi.FieldStorage()
    cgitb.enable()
    sys.stderr = sys.stdout
    print("Content-Type: text/html\n")
    
    main( form ) 
    sys.exit(0)
