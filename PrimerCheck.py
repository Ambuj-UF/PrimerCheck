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


import os
import re
import glob
import argparse
import textwrap
import UserString
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline


parser = argparse.ArgumentParser(prog='PrimerCheck',
                                 version= 'PrimerCheck-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t\t\t Welcome to PrimerCheck-1.0
    \t\t\t\t Use "python PrimerCheck.py -h" for help
    \t\t\t Designed at Kimbal-Braun Lab Group, University of Florida
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))

group = parser.add_mutually_exclusive_group()

group.add_argument('-fq', type=str, default = None, help='Enter fastq file name')
group.add_argument('-fa', type=str, default = None,  help='Enter fasta file name')
parser.add_argument('-frd', type=str, required = True, help='Enter forward primer file name')
parser.add_argument('-rev', type=str, required = True, help='Enter reverse primer file name')
parser.add_argument('-qual', type=str, help='Enter quality data file name')
parser.add_argument('-tag', action='store_true', default=False,
                    help='Use only if intial blastx tagfile output is available in the directory')
parser.add_argument('-tfname', type=str, help='Enter blast tagfile')

args = parser.parse_args()

if args.fa == None and args.fq == None:
    parser.error('fasta or fastq input file required.')

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

if args.tag == True and not args.tfname:
    parser.error('-tfname argument is required in "-tag" mode.')


def primerMatch(record, recordP, ptype):
    
    """Simulation code for mapping primers on the nucleotide sequence and trimming out
        the terminal regions from start to the end of forward primer and from the begening
        of reverse primer to the end of sequence"""
    
    for i, val in enumerate(recordP):
        if ptype == 'forward':
            primerToUse = str(val.seq) if str(val.seq) in str(record.seq) or str(val.seq)[10:-1] in str(record.seq)\
                else str(val.seq)[::-1] if str(val.seq)[::-1] in str(record.seq) or str(val.seq)[::-1][10:-1] in str(record.seq)\
                    else None

        elif ptype == 'reverse':
            primerToUse = str(val.seq.reverse_complement()) if str(val.seq.reverse_complement()) in str(record.seq)\
                or str(val.seq.reverse_complement())[0:15] in str(record.seq) else str(val.seq.reverse_complement())[::-1]\
                    if str(val.seq.reverse_complement())[::-1] in str(record.seq) or str(val.seq.reverse_complement())[::-1][0:15] in str(record.seq)\
                        else None

        if primerToUse == None:
            continue
        else:
            startPosList = [x.start() for x in re.finditer(primerToUse, str(record.seq))]
            if startPosList == [] and ptype == 'reverse':
                startPosList = [x.start() for x in re.finditer(primerToUse[0:15], str(record.seq))]
            if startPosList == [] and ptype == 'forward':
                startPosList = [x.start() for x in re.finditer(primerToUse[10:-1], str(record.seq))]
            endPosList = [x.end() for x in re.finditer(primerToUse, str(record.seq))]
            if endPosList == [] and ptype == 'reverse':
                endPosList = [x.end() for x in re.finditer(primerToUse[0:15], str(record.seq))]
            if endPosList == [] and ptype == 'forward':
                endPosList = [x.end() for x in re.finditer(primerToUse[10:-1], str(record.seq))]
            record.annotations['Primer_' + str(val.id).split('_')[2]] = (val.id)
            return record.annotations, startPosList, endPosList


def blastxOR(BlastInput):
    """Run BLASTX on OR database"""
    if type(BlastInput) is str:
        handleX = open(BlastInput, 'rU'); recordX = list(SeqIO.parse(handleX, 'fasta')); initFlag = False
        tagfile = 'sequenceTag.txt'
    elif type(BlastInput) is list:
        recordX = BlastInput; newRecord = []; initFlag = True; negIDs = []
        tagfile = 'sequenceTagCheck.txt'

    try:
        os.mkdir('OR-Output')
    except OSError:
        pass
                   
    with open(tagfile, 'w') as fp:
        for i, rec in enumerate(recordX):
            stringNuc = ''
            for nuc in rec.seq:
                stringNuc = stringNuc + nuc
        
            blastx_cline = NcbiblastxCommandline(db="ORaaseqs_gallus_taeniopygia.txt",\
                                                evalue=0.0000000001, outfmt=5, out=("OR-Output/Result.%s"%rec.id.split('/')[1]))
            try:
                stdout, stderr = blastx_cline(stdin = stringNuc)
            except:
                continue
                   
            """Parse Blast output"""
            with open('OR-Output/Result.' + str(rec.id).split('/')[1],'r') as xml:
                cFlag = False; eFlag = False; e_val = []; negFlag = False; fcount = 0
                for line in xml:
                    if re.search('No hits found', line) == None:
                        """Check if the sequence belong to OR group"""
                        cFlag = True
                
                    if re.search('<Hsp_evalue>', line) != None:
                        """Extract evalue"""
                        line = line.strip(); line = line.rstrip()
                        line = line.strip('<Hsp_evalue>'); line = line.strip('</')
                        e_val.append(line)
                        eFlag = True
                    
                    if re.search('<Hsp_query-frame>', line) != None and fcount < 1:
                        """Extract frame value"""
                        fcount = fcount + 1
                        line = line.strip(); line = line.rstrip()
                        line = int(line.split('>')[1].split('<')[0])
                        if line < 0:
                            rec.seq = rec.seq.reverse_complement()
                            negFlag = True
            
                if cFlag == True and eFlag == True and negFlag == False:
                    fp.write("%s : OR = True, evalue = %s\n"%(rec.id, e_val[0]))
                elif cFlag == True and eFlag == True and negFlag == True:
                    fp.write("%s : OR = True, evalue = %s, Frame = Negative\n"%(rec.id, e_val[0]))
                else:
                    fp.write("%s : OR = False\n"%rec.id)
                    
                if initFlag == True and cFlag == True and eFlag == True:
                    newRecord.append(rec)
                
                if cFlag == False and initFlag == True and eFlag == False:
                    negIDs.append(rec.id)
    
            os.remove("OR-Output/Result." + str(rec.id).split('/')[1])

    if initFlag == True:
        return newRecord, negIDs


def main():
    if args.fq != None:
        """Executes if Fastq input file supplied"""
        file = args.fq
        with open('sequences.fasta', 'w') as fp, open('RawQual.qual', 'w') as fq:
            handle = open(file, 'rU')
            print("Importing fastq file. This will take a while\n")
            rec = list(SeqIO.parse(handle, 'fastq'))
            print("FastQ file imported\n")
            print("Writing Fasta and Qual data\n")
            SeqIO.write(rec, fp, 'fasta'); SeqIO.write(rec, fq, 'qual')
        
        print("All done. Now importing records\n")
        handle = open('sequences.fasta', 'rU'); records = list(SeqIO.parse(handle, 'fasta'))
        handleQual = open('RawQual.qual', 'rU'); recordsQual = list(SeqIO.parse(handleQual, 'qual'))

    elif args.fa != None:
        """Executes if Fasta input file supplied"""
        handle = open(args.fa, 'rU'); records = list(SeqIO.parse(handle, 'fasta'))
#        handleQual = open(args.qual, 'rU'); # recordsQual = list(SeqIO.parse(handleQual, 'qual'))

    handleF = open(args.frd, 'rU'); recordF = list(SeqIO.parse(handleF, 'fasta'))
    handleR = open(args.rev, 'rU'); recordR = list(SeqIO.parse(handleR, 'fasta'))

    posDict = dict()

    if args.tag == True:
        print("Removing non-OR sequences from the input file\n")
        fopen = open(args.tfname, 'r')
        tagData = fopen.readlines()
        noneIds = [lines.split(' : ')[0] for lines in tagData if 'OR = False' in lines]
        records = [x for x in records if x.id not in noneIds]
        negFrameSeq = [lines.split(' : ')[0] for lines in tagData if 'Frame = Negative' in lines]
        for i, rec in enumerate(records):
            if rec.id in negFrameSeq:
                rec.seq = rec.seq.reverse_complement()
            records[i] = rec

    else:
        print("Running initial Blastx scan...\n")
        records, negIDs = blastxOR(records)
        print("Initial blastx scan completed!\n")

    print("Initiating primer mapping module\n")

    for i, val in enumerate(records):
        fFlag = False; rFlag = False
        try:
            fpAnnotate, startPosListF, endPosListF = primerMatch(val, recordF, 'forward')
        except TypeError:
            fFlag = True
            endPosListF = ['NA']

        try:
            rpAnnotate, startPosListR, endPosListR = primerMatch(val, recordR, 'reverse')
        except TypeError:
            startPosListR = ['NA']
            rFlag = True

        posDict[val.id] = ([endPosListF[0], startPosListR[-1], len(val.seq)])

        if fFlag == False and rFlag == False:
            if endPosListF[0]/len(val.seq) < 0.7:
                records[i].seq = val.seq[endPosListF[0] + 1: startPosListR[-1]]
            else:
                records[i].seq = val.seq[endPosListR[0] + 1: startPosListF[-1]]
        elif fFlag == True and rFlag == False:
            if startPosListR[-1]/len(val.seq) < 0.7:
                records[i].seq = val.seq[endPosListR[0] + 1: -1]
            else:
                records[i].seq = val.seq[0: startPosListR[-1]]
        elif fFlag == False and rFlag == True:
            if endPosListF[0]/len(val.seq) < 0.7:
                records[i].seq = val.seq[endPosListF[0] + 1: -1]
            else:
                records[i].seq = val.seq[0: startPosListF[-1]]
        else:
            pass

    print("Primer mapping done!!\n")


#        for j, recQual in enumerate(recordsQual):
#            if recQual.id == records[i].id:
#                if fFlag == False and rFlag == False:
#                    recQual.letter_annotations["phred_quality"] = recQual.letter_annotations["phred_quality"][endPosListF[0] + 1: startPosListR[-1]]\
#                                                              + [float(x) for x in str('0'*(len(recQual.letter_annotations["phred_quality"])\
#                                                              - len(recQual.letter_annotations["phred_quality"][endPosListF[0] + 1: startPosListR[-1]\
#                                                              ])))]
#                elif fFlag == True and rFlag == False:
#                    recQual.letter_annotations["phred_quality"] = recQual.letter_annotations["phred_quality"][0: startPosListR[-1]]\
#                                                              + [float(x) for x in str('0'*(len(recQual.letter_annotations["phred_quality"])\
#                                                              - len(recQual.letter_annotations["phred_quality"][0: startPosListR[-1]\
#                                                              ])))]
#                elif fFlag == False and rFlag == True:
#                    recQual.letter_annotations["phred_quality"] = recQual.letter_annotations["phred_quality"][endPosListF[0] + 1: -1]\
#                                                              + [float(x) for x in str('0'*(len(recQual.letter_annotations["phred_quality"])\
#                                                              - len(recQual.letter_annotations["phred_quality"][endPosListF[0] + 1: -1\
#                                                              ])))]
#
#            recordsQual[j] = recQual
#
#        if rFlag == False:
#            records[i].annotations.update(dict(rpAnnotate.items()))

    with open('BlastInput.fas', 'w') as fp: #open('ProcessedQual.qual', 'w') as fq:
        SeqIO.write(records, fp, 'fasta')
#        SeqIO.write(recordsQual, fq, 'qual')

    with open('Annotations.txt', 'w') as fp:
        for val in records:
            fp.write('%s: %s\t\t%s\t\t%s\t\t%s\n' %(val.id, val.annotations, posDict[val.id][0], posDict[val.id][1], posDict[val.id][2]))


    print("Initiating second stage blastx search\n")
    blastxOR('BlastInput.fas')
    print("All Done. OR tags are stored in sequenceTag.txt file\n")


                   
if __name__ == "__main__":
    main()















