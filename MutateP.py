################################################################################################################
# Tool for creating primer sequences by substituting nucleotides at variable positions                         #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Braun lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################



import argparse
import textwrap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


parser = argparse.ArgumentParser(prog='MutateP',
                                 version= 'MutateP-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t\t\t Welcome to MutateP-1.0
    \t\t\t\t Use "python PrimerCheck.py -h" for help
    \t\t\t Designed at Kimbal-Braun Lab Group, University of Florida
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))


"""This program takes primer sequence as input and produces all the combinations of primers as per the ambiguous nucleotides present in the primer sequence"""

parser.add_argument('-i', type=str, required = True, help='Enter input primer file name')
parser.add_argument('-o', type=str, required = True, help='Enter output file name')

args = parser.parse_args()

groupDict = {'Y': ['C', 'T'],
    'R': ['A', 'G'],
    'N': ['A', 'C', 'T', 'G'],
    'H': ['A', 'C', 'T'],
    'D': ['A', 'G', 'T'],
    'S': ['C', 'G'],
    'M': ['A', 'C']
    }

def mutate(myList, group):
    retSeqList = [] if group in myList[0] else myList
    if not retSeqList:
        for sequences in myList:
            newSeqList = [sequences.replace(group, nuc, 1) for nuc in groupDict[group] if group in sequences]
            if sequences.count(group) == 0:
                retSeqList = newSeqList
            else:
                for inSeq in newSeqList:
                    for nuc in groupDict[group]:
                        if group in inSeq:
                            newSeqList.append(inSeq.replace(group, nuc, 1))
                for val in newSeqList:
                    if group not in val:
                        retSeqList.append(val)

    else:
        pass

    return retSeqList

def main():
    handle = open(args.i, 'rU')
    records = list(SeqIO.parse(handle, 'fasta'))

    newRecord = []; procRecord = []

    for rec in records:
        sequence = str(rec.seq); myList = [sequence]
        for groups in groupDict.keys():
            myList = mutate(myList, groups)
    
        newRecord.append([SeqRecord(Seq(str(x), IUPACAmbiguousDNA()), id=rec.id + str(i), name=rec.id + str(i),\
                                description=rec.id + str(i)) for i, x in enumerate(myList)])

    for val in newRecord:
        for inval in val:
            procRecord.append(inval)

    with open(args.o, 'w') as fp:
        SeqIO.write(procRecord, fp, 'fasta')


if __name__ == "__main__":
    main()



















