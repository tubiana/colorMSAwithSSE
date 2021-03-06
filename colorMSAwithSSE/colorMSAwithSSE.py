from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio import AlignIO
import Bio
import matplotlib.pyplot as plt
import re
import sys
import numpy as np
from matplotlib import colors
import argparse


class Sequence():

    def __repr__(self, ):
        print(f"ID = {self.id}")
        print(f"sequence (l={len(self.sequence)}):\n{self.sequence}")
        print(f"SSE =\n {self.sse}")
        print(f"alignement (l={len(self.alignment)}) =\n{self.alignment}")
        print(f"SSE aligned =\n{self.sse_aligned}")
        print(f"SSE num =\n{self.sse_num}")
        return ""

    def __init__(self, biopython_ali, structures_folder, alirange):
        self.id = biopython_ali.id
        self.structures_folder = structures_folder

        self.alignment = biopython_ali.seq
        self.sequence = self.alignment.ungap()
        self.pdbname = None
        self.sse = None
        self.sse_ali = None
        self.non_empty_position = []
        self.strip_start = 0
        self.alignment_sequence_matching = {}

        self.define_position_in_ali()
        self.define_pdbname()
        self.pdbpath = f"{self.structures_folder}/{self.pdbname}.pdb"

        self.add_pdb()
        self.add_sse()
        self.set_sse_as_number()
        if alirange != 'all':
            self.strip(alirange)



    def define_position_in_ali(self):
        gap_char = ["_", "-"]

        sequence_position = 0
        for i in range(len(self.alignment)):
            letter = self.alignment[i]
            self.alignment_sequence_matching[i]=sequence_position
            if letter not in gap_char:
                self.non_empty_position.append(i)
                sequence_position+=1

    def define_pdbname(self):
        regexCATH = re.compile("([0-9][\w]{3}[A-Z][0-9]{2})")
        regexPDB = re.compile("([0-9][\w]{3})")
        regexPDBfile = re.compile("(\w*)\.pdb")
        matchCATH = regexCATH.match(self.id)
        matchPDB = regexPDB.match(self.id)
        matchPDBfile = regexPDBfile.match(self.id)
        if not matchCATH and not matchPDB and not matchPDBfile:
            print(f"NO PDB ID FOUND IN HEADER ID FOR SEQUENCE {self.id}. PLEASE CHECH AND CORRECT")
        if matchCATH:
            pdbname = matchCATH.group(1)
        elif matchPDB:
            pdbname = matchPDB.group(1)
        elif matchPDBfile:
            pdbname = matchPDBfile.group(1)
        self.pdbname = pdbname

    def add_pdb(self):
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure(self.pdbname, self.pdbpath)[0]

    def add_sse(self, ):

        def simplify_SS(ss):
            """
            Secondary Structure simplification from DSSP.
            All helix types are just helix
            All sheet type are just sheet
            turn, bend and coil are 'C' (coil)
            Args:
                ss (str): secondary structure letter
            Returns:
                str: simplified secondary structure letter.
            """
            if ss in ["H", "G", "I"]:
                return "H"
            elif ss in ["E", "B"]:
                return "E"
            elif ss in ["T", "S", "-"]:
                return "C"

        dssp = DSSP(self.structure, self.pdbpath, dssp="mkdssp")
        self.sse_full = "".join([dssp[x][2] for x in dssp.keys()])
        self.sse = "".join([simplify_SS(dssp[x][2]) for x in dssp.keys()])

        lenseq = len(self.non_empty_position)
        lendssp = len(self.sse)

        if lenseq != lendssp:
            print(f"ERROR, secondary structures {lenseq} and sequence {lendssp} don't have the same size.")
            print(f" File: {self.pdbpath}")
            print(f"This could be due to: ")
            print(f" - mistmatch between the sequence and secondary structures")
            print(f" - DSSP that ignored some CTER or NTER residues (try DSSP on this file to check)")
            print(f" - Duplicated residues in your sequence.")
            sys.exit(1)

        sse_aligned = ["-"] * len(self.alignment)
        for i in range(len(self.non_empty_position)):
            sse_aligned[self.non_empty_position[i]] = self.sse[i]
        self.sse_aligned = "".join(sse_aligned)

    def set_sse_as_number(self):
        conversion = {"-": 0,
                      "C": 1,
                      "H": 2,
                      "E": 3}
        self.sse_num = [conversion[x] for x in self.sse_aligned]

    def strip(self, alirange):
        startali = int(alirange.split(':')[0]) -1
        self.strip_start = startali
        endali = int(alirange.split(':')[1])
        if startali > endali:
            print("ERROR while striping sequence, starting position is higher than ending position")
            sys.exit(1)
        startseq = self.alignment_sequence_matching[startali]
        endseq = self.alignment_sequence_matching[endali]

        self.alignment = self.alignment[startali:endali]
        self.sse_aligned = self.sse_aligned[startali:endali]
        self.sse_num = self.sse_num[startali:endali]

        self.sequence = self.sequence[startseq:endseq]
        self.sse_full = self.sse_full[startseq:endseq]
        self.sse = self.sse[startseq:endseq]





def parseArg():
    """
    Parse the arguments
    return:
        args (dict): dictionnary of arguments

    """
    arguments = argparse.ArgumentParser(description="This program is made to generate a picture (SVG/PNG/PDF) of a multiple"
                                                    " sequence alignement with secondary structure element coloration",
                                        usage="python colorMSAwithSSE.py -f file.fasta -t structuresfolder -o alignment.png")
    try:
        argcomplete.autocomplete(arguments)
    except:
        pass
    arguments.add_argument('-f', "--file", help="alignment file (FASTA only for now)", required=True)
    arguments.add_argument('-s', '--structures', help="structures FOLDER", required=True)
    arguments.add_argument('-w', '--linewidth', help="line width. Default = 50, choose 0 if you want the picture one 1 line", default=50, type=int)
    arguments.add_argument('-i', '--interspace', help="interspace between two lines", default=2, type=int)
    arguments.add_argument('-fs', '--fontstyle', help="fontstyle (serif/sans). Default: serif.", default="serif")
    arguments.add_argument('-o', '--output', help="output file", default="out.png")
    arguments.add_argument('-sort', '--sort', help="Sort indexes (Y/N) (A->Z)", default="Y")
    arguments.add_argument('-r', '--range', help="output range of the alignment. Example '30:70'. default='all'", default="all")
    arguments.add_argument('-p', '--tickPosition', help="Put the position in the alignment every {tickPosition} ", default=10, type=int)
    arguments.add_argument('-c','--color', help="Color for the secondary structure in this order 'GAP,COIL,HELICES,BSHEET'. "
                           "default is 'white, bisque, red, yellow'. You can use color name (See https://matplotlib.org/stable/gallery/color/named_colors.html) or HEX code", default="white,bisque,red,yellow")
    arguments.add_argument('-a', '--alpha', help="Transparency of the color background", default=0.6, type=float)




    args = vars(arguments.parse_args())
    if args["sort"] == "Y":
        args["sort"]=True
    return (args)


def read_alignment(path:str):
    """
    Read an alignment from a filepath and define it's type.
    :param path: Path of the alignment file
    :return: Biopython Alignment object
    """
    alitype = {'fa':'fasta',
               'faa':'fasta',
               'fasta':'fasta',
               'clw':'clustal',
               'aln':'clustal',
               'clustal':'clustal',
               'phy':'phylip',
               'phylip':'phylip',
               }

    extension = path.split('.')[-1]
    if extension not in alitype:
        print("alignment Extension not in the supported list for now")
        print("Please convert your file into fasta forma (with 'fa','faa' or 'fasta' extension)")
        print("or open a ticket on this github repository https://github.com/tubiana/colorMSAwithSSE")
        sys.exit(1)

    return AlignIO.read(path, alitype[extension])
    
    return 


def create_ali_list(alignment:Bio.Align.MultipleSeqAlignment, structures_folder:str, sort=True, alirange='all'):
    """
    Generate a list of alignment objects from the MSA file
    :param alignment: Biopython alignment Object
    :param structures_folder: Folder location where the structures are
    :param sort: Sort labels (pdbid).
    :param alirange: alignment range.
    :return seqlist: list of Sequence object
    """
    seqlist = []
    for ali in alignment:
        seqlist.append(Sequence(ali, structures_folder, alirange))
        
    #Order alphabetically the list
    if sort==True:
        seqlist.sort(key=lambda x: x.pdbname)
    return seqlist


def chunk(data:list, N:int):
    """
    Chunk a list into a list of list of size N
    :param data: list to chunk
    :param N: Chunk size
    :return: Chucked list (list of list)
    """
    array = []
    for i in range(0, len(data), N):
        line = data[i:i + N]
        array.append(line)
    return array


def prepare_data(seqlist:list, chunksize:int, interspace:int, tickPosition=10, strip_position=(0,0)):
    """
    Prepare input matrix and labels for the plot.
    :param seqlist: list of sequence objects
    :param chunksize: chunk size
    :param interspace: space between 2 group of sequences
    :param tickPosition: Label the position every X
    :param strip_position: if the sequence is cut (or 'stripped'), this value will adjust the position on the graph.
    :return datamatrix: Matrix with all SSE numbers
    :return labels_lines: Labels lines
    :return data_annot: Matrix with all annotations (sequences)
    :return data_positions: Matrix with all position numbers
    """


    #1. Instantiate the dictionnaries that will contain the SSE and the sequence
    chunked_sse = {}
    chunked_seq = {}


    #1. Prepare the position list (tick every X position
    ALISIZE = len(seqlist[0].alignment)
    chunked_position = [""] * ALISIZE
    for i in range(0, ALISIZE, tickPosition):
        chunked_position[i] = i + 1 + strip_position[0]
    if strip_position[1] != 0:
        chunked_position[-1] = strip_position[1]
    else:
        chunked_position[-1] = ALISIZE

    #2. Chunk lists
    for i in range(len(seqlist)):
        if chunksize > 0:
            chunked_sse[seqlist[i].pdbname] = chunk(seqlist[i].sse_num, chunksize)
            chunked_seq[seqlist[i].pdbname] = chunk(list(seqlist[i].alignment), chunksize)

        else:
            chunked_sse[seqlist[i].pdbname] = [seqlist[i].sse_num]
            chunked_seq[seqlist[i].pdbname] = [list(seqlist[i].alignment)]

    if chunksize > 0:
        chunked_position = chunk(chunked_position, chunksize)
        WIDTH = chunksize
    else:
        chunked_position = [chunked_position]
        WIDTH = ALISIZE

    #3. Instiation of objects that will contains the data and labels.
    labels = list(chunked_sse.keys())
    labels_lines = []
    data_matrix = []
    data_annot = []
    data_positions = []

    #3.1 fill the first line as empty lines because it will contains the positions
    data_matrix.append([0] * WIDTH)
    data_annot.append([""] * WIDTH)

    #4. Define the number of lines and get the number of sequences
    NSEQ = len(labels)
    NLINES = len(list(chunked_sse.values())[0])

    #5. Fill the matrices with numbers for SSE and letters for labels.
    # For every lines
    for i in range(NLINES):
        #The first line will always be the position
        line_position = chunked_position[i]
        #And then we have to be carefull that the last one is fill with 0 to have the same lines width everywere.
        if len(line_position) < WIDTH:
            data_positions.append(line_position + [""] * (WIDTH - len(line_position)))
            labels_lines.append("Position")
        else:
            data_positions.append(line_position)
            labels_lines.append("Position")

        #Then for every sequence
        for j in range(NSEQ):
            #We get the labels/SSE numbers
            seq = labels[j]
            labels_lines.append(seq)
            line = chunked_sse[seq][i]
            lineseq = chunked_seq[seq][i]

            curr_linewidth = len(line) #we get the current linewidth

            # The last block is usually smaller than other chunked block so we need to fill it with empty values
            # because imshow only takes matrices
            if curr_linewidth < WIDTH:
                line = line + [0] * (WIDTH - curr_linewidth)
                lineseq = lineseq + [""] * (WIDTH - curr_linewidth)

            #Add the lines into the matrix (list of list)
            data_matrix.append(line)
            data_annot.append(lineseq)
            data_positions.append([""] * WIDTH)

        #If it's not the last line. Add the "Interspace"
        if not i == NLINES - 1:
            for _ in range(interspace):
                data_matrix.append([0] * WIDTH)
                data_annot.append([""] * WIDTH)
            for _ in range(interspace - 1):
                data_positions.append([""] * WIDTH)
                labels_lines.append("")

    #Convert the SSE elements numbers to a numpy matrix for matplotlib
    data_matrix = np.asarray(data_matrix)

    return (data_matrix,
            labels_lines,
            data_annot,
            data_positions,)


def generate_graph(data_matrix:np.array, labels_lines:list, data_annot:list,
                   data_positions:list, output:str, color:list, alpha:float, fontstyle:str,
                   alirange:str):
    """
    Generate the figure
    :param datamatrix: Matrix with all SSE numbers
    :param labels_lines: Labels lines
    :param data_annot: Matrix with all annotations (sequences)
    :param data_positions: Matrix with all position numbers
    :param output: output filename
    :param colors: colors list in this order Gap, Coil, Helices, BSheet
    :param alpha: Transparency
    :param fontstyle: font style (serif or sans-serif)
    :param alirange: alignment range in this syntax "XX-YY"
    :return: None
    """

    if fontstyle == "serif":
        plt.rcParams["font.family"] = "DejaVu Serif"
    else:
        plt.rcParams["font.family"] = "DejaVu Sans"

    WIDTH = len(data_matrix[0])
    NLINES = len(data_matrix)

    fig, ax = plt.subplots(figsize=(WIDTH/3, WIDTH/3), facecolor='w')


    cmap = colors.ListedColormap(color)
    boundaries = list(range(len(color) + 1))
    norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

    plt.imshow(data_matrix,
               interpolation=None,
               cmap=cmap,
               norm=norm,
               alpha=alpha,
               )

    #Add labels
    for i in range(NLINES):
        for j in range(WIDTH):
            label = data_annot[i][j]
            label_position = data_positions[i][j]
            if not label == "":
                ax.text(j, i, label, color='black', ha='center', va='center', size="12")
            if not label_position == "":
                ax.text(j, i, label_position, color='black', ha='center', va='center', size="8")
                ax.text(j, i+0.5, "|", color='black', ha='center', va='center', size="6")

    # Set Xticks
    ytick_position = []
    ytick_labels = []
    for i in range(len(labels_lines)):
        if labels_lines[i] != "":
            ytick_position.append(i)
            ytick_labels.append(labels_lines[i])
    ax.set_yticks(ytick_position)
    ax.set_yticklabels(ytick_labels)
    ax.set_xticks([])
    ax.set_xticklabels([])

    plt.savefig(output,
                dpi=300, bbox_inches='tight',
                pad_inches=0)
    plt.close()


def get_alirange(alirange):
    ali_regex = re.compile("(all)|((\d+):(\d+))")
    match = ali_regex.math(alirange)
    if not match:
        print("Error on alignment range. Please use 'all' for all alignment")
        print("or 'BEGIN:END' to strip one part of the sequence")
        sys.exit(1)
    return alirange

def engine():
    """
    Main function
    :return: None
    """

    args = parseArg()

    alignmentfile = args["file"]
    structures_folder = args["structures"]

    chunksize = args["linewidth"]
    interspace = args["interspace"]
    output = args["output"]
    tickPosition = args["tickPosition"]
    fontstyle=args["fontstyle"]
    sort = args["sort"]
    alirange=args["range"]
    color=[x.strip() for x in args["color"].split(',')]
    alpha = args["alpha"]

    alignment = read_alignment(alignmentfile)
    seqlist = create_ali_list(alignment, structures_folder, sort, alirange)

    if alirange != 'all':
        strip_position = (int(alirange.split(':')[0]) -1,
                          int(alirange.split(':')[1]) )
    else:
        strip_position = (0,0)

    (data_matrix,
     labels_lines,
     data_annot,
     data_positions) = prepare_data(seqlist, chunksize, interspace, tickPosition,strip_position)

    print("Generating Graphs.. Please wait...")
    generate_graph(data_matrix=data_matrix,
                   labels_lines=labels_lines,
                   data_annot=data_annot,
                   data_positions=data_positions,
                   output=output,
                   color=color,
                   alpha=alpha,
                   fontstyle=fontstyle,
                   alirange=alirange)

    print(f"> Done. Check {output}")





if __name__ == "__main__":
    engine()


