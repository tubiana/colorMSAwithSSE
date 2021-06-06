# colorMSAwithSSE
Generate an image from a multiple sequence alignment (MSA) colored by the secondary structures elements (SSE)

## Usage
**Please note that you need a strict match between the sequences in the fasta file AND your structrure files**. It will not work if there is a mismatch.

1. You will need a fasta where the header countain the PDB filename
2. And all the structures used

### formats
For now only those file extension and associated formats are supported:
- `.fa` -> fasta
- `.faa` -> fasta
- `.fasta` -> fasta
- `.clw` -> clustal
- `.aln` -> clustal
- `.clustal` -> clustal
- `.phy` -> phylip
- `.phylip` -> phylip

### Coloring
You can change the color of the Gaps, Coil, ɑ-helix and β-sheet with the argument `-c`.  
Use the HEX code or the matplotlib color name here (https://matplotlib.org/stable/gallery/color/named_colors.html).  
The color order should follow this : GAP, COIL, HELICES, BSHEET.  
Example: `-c "white, bisque, red, yellow"`

### Font
You can change the font style with the parameter `-fs`. It can be with Serif (`serif`) or sans-serif (`sans`)


### Example
You can run an example within the test folder.
  
`python ../colorMSAwithSSE/colorMSAwithSSE.py -f SH2.fasta -s structures/ -o SH2.png -w 60 -c 'white, bisque, red, yellow' -a 0.6 -ft sans`

**Output example**
![outputexample](test/SH2.png) 


## help

```
usage: python colorMSAwithSSE.py -f file.fasta -t structuresfolder -o alignment.png

This program is made to generate a picture (SVG/PNG/PDF) of a multiple
sequence alignement with secondary structure element coloration

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  alignment file (FASTA only for now)
  -s STRUCTURES, --structures STRUCTURES
                        structures FOLDER
  -w LINEWIDTH, --linewidth LINEWIDTH
                        line width. Default = 50, choose 0 if you want the
                        picture one 1 line
  -i INTERSPACE, --interspace INTERSPACE
                        interspace between two lines
  -fs FONTSTYLE, --fontstyle FONTSTYLE
                        fontstyle (serif/sans). Default: serif.
  -o OUTPUT, --output OUTPUT
                        output file
  -sort SORT, --sort SORT
                        Sort indexes (Y/N) (A->Z)
  -p TICKPOSITION, --tickPosition TICKPOSITION
                        Put the position in the alignment every {tickPosition}
  -c COLOR, --color COLOR
                        Color for the secondary structure in this order
                        'GAP,COIL,HELICES,BSHEET'. default is 'white, bisque,
                        red, yellow'. You can use color name (See https://matp
                        lotlib.org/stable/gallery/color/named_colors.html) or
                        HEX code
  -a ALPHA, --alpha ALPHA
                        Transparency of the color background
```
## Dependancies:
- Biopython
- matplotlib
- DSSP
- argparse

You can the dependancies with `conda install -c conda-forge -c salilab dssp biopython matplotlib`


## Support
Please ask your question/feedback through Github :-)
