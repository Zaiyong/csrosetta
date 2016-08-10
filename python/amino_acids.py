import string
from operator import add

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u'}

molecular_weight={'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}

def long2short(name):
    return longer_names[name.upper()]


short_to_long = {}
for rsd in longer_names.keys():short_to_long[longer_names[rsd]] = rsd

def short2long(name):
    return short_to_long[name.upper()]

hpcg = {'A': 'H', 'V': 'H', 'I': 'H', 'L': 'H', 'F': 'H', 'P': 'H', 'M': 'H',
        'G': 'G',
        'D': 'C', 'E': 'C', 'K': 'C', 'R': 'C',
        'S': 'P', 'T': 'P', 'Y': 'P', 'C': 'P', 'N': 'P', 'Q': 'P', 'H': 'P',
        'W': 'P'}

size = {'G': 66, 'A': 92, 'V': 142, 'L': 168, 'I': 169,
        'S': 99, 'T': 122, 'C': 106, 'P': 129, 'F': 203,
        'Y': 204, 'W': 238, 'H': 167, 'D': 125, 'N': 135,
        'E': 155, 'Q': 161, 'M': 171, 'K': 171, 'R': 202,
        'X': 10000}

frequency = {'G': .075, 'A': .09, 'V': .069, 'L': .075, 'I': .046,
               'S': .071, 'T': .060, 'C': .028, 'P': .046, 'F': .035,
               'Y': .035, 'W': .011, 'H': .021, 'D': .055, 'N': .044,
               'E': .062, 'Q': .039, 'M': .017, 'K': .07, 'R': .047}

HP = {'I': 0.73, 'F': 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37,
      'M': 0.26, 'A': 0.25, 'G': 0.16, 'C': 0.04, 'Y': 0.02,
      'P': -0.07, 'T': -0.18, 'S': -0.26, 'H': -0.40, 'E': -0.62,
      'N': -0.64, 'Q': -0.69, 'D': -0.72, 'K': -1.10, 'R': -1.76}

HP['X'] = reduce(add,HP.values())/20
GES = {'F': -3.7, 'M': -3.4, 'I': -3.1, 'L': -2.8, 'V': -2.6,
       'C': -2.0, 'W': -1.9, 'A': -1.6, 'T': -1.2, 'G': -1.0,
       'S': -0.6, 'P': 0.2,  'Y': 0.7,  'H': 3.0,  'Q': 4.1,
       'N': 4.8,  'E': 8.2,  'K': 8.8,  'D': 9.2,  'R': 12.3}

GES['X'] = reduce(add,GES.values())/20


#from http://astral.stanford.edu/scopseq-1.55/release-notes-1.55.txt
extra_longer_names={'TYY': 'Y', 'ARM': 'R', 'ARG': 'R', 'HTR': 'W', 'C5C': 'C', 'CY1': 'C', 'DPN': 'F', 'BHD': 'D', 'SCS': 'C', 'TPQ': 'A', 'MEN': 'N', 'SEP': 'S', 'KCX': 'K', 'TYB': 'Y', 'NMC': 'G', 'CSS': 'C', 'CSP': 'C', 'SET': 'S', 'EFC': 'C', 'CSW': 'C', 'TRP': 'W', 'BMT': 'T', 'CSX': 'C', 'SHC': 'C', 'NLE': 'L', 'OCS': 'C', 'DSN': 'S', 'DHA': 'A', 'TPL': 'W', 'HPQ': 'F', 'CSD': 'A', 'MSE': 'M', 'MLE': 'L', 'DLE': 'L', 'CSO': 'C', 'CY3': 'C', 'NLN': 'L', 'CYM': 'C', 'BUG': 'L', 'SAR': 'G', 'PEC': 'C', 'BUC': 'C', 'DAS': 'D', 'DAR': 'R', 'CYG': 'C', 'HIC': 'H', 'HMR': 'R', 'DAH': 'F', 'LLP': 'K', 'DAL': 'A', 'CEA': 'C', 'MPQ': 'G', 'HIS': 'H', 'CYQ': 'C', 'NLP': 'L', 'CYS': 'C', 'AIB': 'A', 'MSA': 'G', 'SEL': 'S', 'PCA': 'E', 'DIV': 'V', 'DTR': 'W', 'FLA': 'A', 'IYR': 'Y', 'DNP': 'A', '2AS': 'D', 'DTY': 'Y', 'HYP': 'P', 'ASX': 'B', 'SOC': 'C', 'CME': 'C', 'PRR': 'A', 'DPR': 'P', 'DTH': 'T', 'DIL': 'I', 'TYS': 'Y', 'TYR': 'Y', 'CXM': 'M', 'GL3': 'G', 'STY': 'Y', 'ALY': 'K', 'LEU': 'L', 'TRG': 'K', 'MET': 'M', 'MAA': 'A', 'ALA': 'A', 'PHE': 'F', 'SHR': 'K', 'HIP': 'H', 'PHL': 'F', 'ALM': 'A', 'PHI': 'F', 'ALO': 'T', 'MVA': 'V', 'OMT': 'M', 'FME': 'M', 'HAR': 'R', 'DLY': 'K', '5HP': 'E', 'BNN': 'A', 'PTR': 'Y', 'SCY': 'C', 'MHS': 'H', 'SVA': 'S', 'TRO': 'W', 'THR': 'T', 'AYA': 'A', 'HAC': 'A', 'GMA': 'E', 'CLE': 'L', 'TPO': 'T', 'CHG': 'A', 'VAL': 'V', 'LYZ': 'K', 'GSC': 'G', 'LTR': 'W', 'ASL': 'D', 'OAS': 'S', 'ASN': 'N', 'CGU': 'E', 'LYM': 'K', 'PR3': 'C', 'GGL': 'E', 'ASK': 'D', 'DCY': 'C', 'SCH': 'C', 'ASA': 'D', 'ASB': 'D', 'DGL': 'E', 'DGN': 'Q', 'IIL': 'I', 'PRO': 'P', 'SAC': 'S', 'LYS': 'K', 'SER': 'S', 'ASP': 'D', 'ASQ': 'D', 'NEM': 'H', 'BCS': 'C', 'NEP': 'H', 'GLU': 'E', 'SMC': 'C', 'DSP': 'D', 'MIS': 'S', 'TYQ': 'Y', 'TIH': 'A', 'C6C': 'C', 'GLZ': 'G', 'GLY': 'G', '3AH': 'H', 'ACL': 'R', 'CCS': 'C', 'DVA': 'V', 'LLY': 'K', 'GLX': 'Z', 'DHI': 'H', 'GLN': 'Q', 'ILE': 'I', 'AGM': 'R', 'PAQ': 'Y'}

##load the G-X-G surface areas
## surface_area = {}
## data = open('residue_surface_areas.dat','r')
## line = data.readline()
## line = data.readline()
## for i in range(20):
##     l = string.split(line)
##     assert l[0] in amino_acids
##     surface_area[l[0]] = float(l[1])
##     line = data.readline()
## data.close()

#surface area from Rost and Sander Proteins 20,216-226
SA = {'A':106, 'C':135, 'D':163, 'E':194, 'F':197,
      'G': 84, 'H':184, 'I':169, 'K':205, 'L':164,
      'M':188, 'N':157, 'P':136, 'Q':198, 'R':248,
      'S':130, 'T':142, 'V':142, 'W':227, 'Y':222}

# number of chemical shift that should have for each amino acid, here O and C and S are not included.

NO_resonances={'ALA': 6, 'ARG': 22, 'ASN': 11, 'ASP': 9, 'CYS': 8,
               'GLN': 14, 'GLU': 12, 'GLY': 5, 'HIS': 16, 'ILE': 13,
               'LEU': 13, 'LYS': 18, 'MET': 12, 'PHE': 18, 'PRO': 13,
               'SER': 8, 'THR': 9, 'TRP': 22, 'TYR': 18, 'VAL': 10}

NO_bb_resonances={'ALA': 5, 'ARG': 5, 'ASN': 5, 'ASP': 5, 'CYS': 5,
               'GLN': 5, 'GLU': 5, 'GLY': 5, 'HIS': 5, 'ILE': 5,
               'LEU': 5, 'LYS': 5, 'MET': 5, 'PHE': 5, 'PRO': 4,
               'SER': 5, 'THR': 5, 'TRP': 5, 'TYR': 5, 'VAL': 5}
