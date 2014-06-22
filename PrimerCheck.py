################################################################################################################
# Tool for mapping primers and annotating fasta and fastq sequences                                            #
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

"""Code requires local protein sequence database/file to perform Blastx operation. 
This executes for OR database which I obtained from Johny Wright, PhD student, 
Biology Department, University of Florida. A slight change in database name supplied 
in line 147 will do well for others"""

try:
    from Bio.Blast.Applications import NcbiblastxCommandline
    from Bio import SeqIO
except: ImportError, e:
    sys.exit("BioPython not found on your system. Program Terminated")


import os
import re
import sys
import argparse
import textwrap


parser = argparse.ArgumentParser(prog='PrimerCheck',
                                 version= 'PrimerCheck-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t\t\t Welcome to PrimerCheck-1.0
    \t\t\t\t Use "python PrimerCheck.py -h" for help
    \t\t\t Designed at Kimball-Braun Lab Group, University of Florida
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))

group = parser.add_mutually_exclusive_group()

group.add_argument('-fq', type=str, default = None, help='Enter fastq file name')
group.add_argument('-fa', type=str, default = None,  help='Enter fasta file name')
parser.add_argument('-frd', type=str, help='Enter forward primer file name')
parser.add_argument('-rev', type=str, help='Enter reverse primer file name')
parser.add_argument('-qual', type=str, help='Enter quality data file name')

args = parser.parse_args()

if args.fa == True and not args.frd:
    parser.error('-frd argument is required in "-fa" mode.')

if args.fa == True and not args.rev:
    parser.error('-rev argument is required in "-fa" mode.')

if args.fa == True and not args.qual:
    parser.error('-qual argument is required in "-fa" mode.')

if args.fq == True and not args.frd:
    parser.error('-frd argument is required in "-fa" mode.')

if args.fq == True and not args.qual:
    parser.error('-rev argument is required in "-fq" mode.')


def oneNucChange(primer):
    """Generates list of primer with one nucleotide substitution"""
    for i, nuc in enumerate(str(primer)):
        lprime = list(primer)
        for val in ['A', 'C', 'G', 'T']:
            myList.append(''.join(lprime[i] = val))
    return myList

def twoNucChange(primer):
    """Generates list of primer with two nucleotide substitution"""
    primers = oneNucChange(primer):
    for primerData in primers:
        for i, nuc in enumerate(str(primerData)):
            lprime = list(primerData)
            for val in ['A', 'C', 'G', 'T']:
                myList.append(''.join(lprime[i] = val))
    return myList

def threeNucChange(primer):
    """Generates list of primer with three nucleotide substitution"""
    primers = twoNucChange(primer):
    for primerData in primers:
        for i, nuc in enumerate(str(primerData)):
            lprime = list(primerData)
            for val in ['A', 'C', 'G', 'T']:
                myList.append(''.join(lprime[i] = val))
    return myList

def prodPrimer(primerList, record):
    """Checks if primer matches with the input nucleotide sequence"""
    for primers in primerList:
        if primers in str(record.seq) == True:
            return primers, True
    return None, False

def primerMatch(record, recordP):
    """Simulation code for mapping primers on the nucleotide sequence and trimming out
        the terminal regions from start to the end of forward primer and from the begening
        of reverse primer to the end of sequence"""
    oneFlag = False
    twoFlag = False
    for i, val in enumerate(recordP):
        primerToUse = str(val.seq) if str(val.seq) in str(record.seq) == True else None
        if primerToUse == None:
            primerToUse, oneFlag = prodPrimer(oneNucChange(val.seq), record)
        if oneFlag == Flase:
            primerToUse, twoFlag = prodPrimer(twoNucChange(val.seq), record)
        if twoFlag == False:
            primerToUse = prodPrimer(threeNucChange(val.seq), record)[0]

        continue if primerToUse == None else pass
        startPosList = [x.start() for x in re.finditer(primerToUse, str(record.seq))]
        endPosList = [x.end() for x in re.finditer(primerToUse, str(record.seq))]
        record.annotations['Primer_' + str(val.id)] = (val.id)
        return record.annotations, startPosList, endPosList

def blastxOR(BlastInput):
    """Run BLASTX on OR database"""
    handleX = open(BlastInput, 'rU')
    recordX = list(SeqIO.parse(handleX, 'fasta'))
    try:
        os.mkdir('OR-Output')
    except OSError:
        pass
                   
    with open('sequenceTag.txt', 'w') as fp:
        for i, rec in enumerate(recordX):
            blastx_cline = NcbiblastxCommandline(query=rec, db="ORaaseqs_gallus_taeniopygia.txt",\
                                                evalue=e-10, outfmt=5, out=("OR-Output/Result.%s"%rec.id))
            stdout, stderr = blastx_cline()
                   
            """Parse Blast output"""
            with open('OR-Output/Result.' + str(rec.id),'r') as xml:
                Flag = False; eFlag = False; e_val = []
                for line in xml:
                    if re.search('No hits found', line) != None:
                        """Check if the sequence belong to OR group"""
                        Flag = True
                    if re.search('<Hsp_evalue>', line) != None:
                        """Extract evalue"""
                        line = line.strip(); line = line.rstrip();
                        line = line.strip('<Hsp_evalue>'); line = line.strip('</')
                        e_val.append(line)
                        eFlag = True
                   
                if Flag == True and eFlag == True:
                    fp.write("%s : OR = True, evalue = %s"%(rec.id, e_val[0]))
                else:
                    fp.write("%s : OR = False"%rec.id)
                   

def main():
    if args.fq != None:
    """Executes if Fastq input file supplied"""
        file = args.fq
        with open('sequences.fasta', 'w') as fp, open('RawQual.qual', 'w') as fq:
            handle = open(file, 'rU')
            print("Importing fastq file. This will take a while\n")
            rec = list(SeqIO.parse(handle, 'fastq-illumina'))
            SeqIO.write(rec, fp, 'fasta'); SeqIO.write(rec, fr, 'qual')
            
        handle = open('sequences.fasta', 'rU'); records = list(SeqIO.parse(handle, 'fasta'))
        handleQual = open('RawQual.qual', 'rU'); recordsQual = list(SeqIO.parse(handleQual, 'qual'))

    elif args.fa != None:
        """Executes if Fasta input file supplied"""
        handle = open(args.fa, 'rU'); records = list(SeqIO.parse(handle, 'fasta'))
        handleQual = open(args.qual, 'rU'); recordsQual = list(SeqIO.parse(handleQual, 'qual'))

    handleF = open(args.frd, 'rU'); recordF = list(SeqIO.parse(handleF, 'fasta'))
    handleR = open(args.rev, 'rU'); recordR = list(SeqIO.parse(handleR, 'fasta'))

    for i, val in enumerate(records):
        try:
            fpAnnotate, startPosListF, endPosListF = primerMatch(records, recordF)
        except TypeError:
            print("Forward primer match not found for %s" %val.id)
            continue
        try:
            rpAnnotate, startPosListR, endPosListR = primerMatch(records, recordR)
        except TypeError:
            print("Reverse primer match not found for %s" %val.id)
            continue

        records[i].seq = val.seq[endPosListF[0] + 1: startPosListF[-1]]
        for j, recQual in enumerate(recordsQual):
            if recQual.id == records[i].id:
                recQual.seq = recQual.seq[endPosListF[0] + 1: startPosListF[-1]]
            recordsQual[j] = recQual
        records[i].annotations.update(dict(fpAnnotate.items() + rpAnnotate.items()))

    with open('BlastInput.fas', 'w') as fp, open('ProcessedQual.qual', 'w') as fq:
        SeqIO.write(records, fp, 'fasta')
        SeqIO.write(recordsQual, fq, 'qual')

    with open('Annotations.txt', 'w') as fp:
        for val in records:
            fp.write('%s\n' %val.annotations)

    blastxOR('BlastInput.fas')


                   
if __name__ == "__main__":
    main()

