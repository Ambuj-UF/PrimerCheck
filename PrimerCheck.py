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
import sys
import glob
import shutil
import argparse
import textwrap
import collections
import UserString
import pylab as pl
import numpy as np
from Bio import SeqIO
from collections import defaultdict
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
parser.add_argument('-tag', action='store_true', default=False,
                    help='Use only if intial blastx tagfile output is available in the directory')
parser.add_argument('-tfname', type=str, help='Enter blast tagfile')
parser.add_argument('-o', type=str, help='Enter output folder name')

args = parser.parse_args()

if args.fa == None and args.fq == None:
    parser.error('fasta or fastq input file required.')

if args.fa == True and not args.frd:
    parser.error('-frd argument is required in "-fa" mode.')

if args.fa == True and not args.rev:
    parser.error('-rev argument is required in "-fa" mode.')

if args.fq == True and not args.frd:
    parser.error('-frd argument is required in "-fa" mode.')

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
    
    seqList1 = []; seqList2 = []; seqList3 = []

    toolbar_width = len(recordX)
    with open(tagfile, 'w') as fp:
        for i, rec in enumerate(recordX):
            p = str((float(i)/toolbar_width)*100)[:4]
            sys.stdout.write("\r%s%%" %p)
            sys.stdout.flush()
            
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



def processOut(outputFolder):
    recordsF1= list()
    recordsF2= list()
    recordsF3= list()
    recordsQ1= list()
    recordsQ2= list()
    recordsQ3= list()
    dataDict = dict()
    dataDictAnnot = dict()
    
    handleFasta = open('Sequences.fasta', 'rU')
    recordsFasta = list(SeqIO.parse(handleFasta, 'fasta'))
    
    handleQual = open('RawQual.qual', 'rU')
    recordsQual = list(SeqIO.parse(handleQual, 'qual'))
    
    fileData = open('sequenceTagCheck.txt', 'r').readlines()
    
    anFile = open('Annotations.txt', 'r').readlines()
    
    for lines in anFile:
        dataDictAnnot[lines.split(' ')[0][0:-1]] = {'Begening': lines.split('\t\t')[1], 'End': lines.split('\t\t')[2], 'length': lines.split('\t\t')[3]}
    
    def seqTrim(records):
        for i, rec in enumerate(records):
            try:
                if dataDictAnnot[rec.id]['Begening'] == 'NA':
                    dataDictAnnot[rec.id]['Begening'] = 0
                if dataDictAnnot[rec.id]['End'] == 'NA':
                    dataDictAnnot[rec.id]['End'] = dataDictAnnot[rec.id]['length']
                
                records[i].seq = rec.seq[int(dataDictAnnot[rec.id]['Begening']): int(dataDictAnnot[rec.id]['End'])]
            
            except KeyError:
                continue
        
        return records
    
    
    def qualEdit(records):
        newRecords = list()
        for i, rec in enumerate(records):
            try:
                if dataDictAnnot[rec.id]['Begening'] == 'NA':
                    dataDictAnnot[rec.id]['Begening'] = 0
                if dataDictAnnot[rec.id]['End'] == 'NA':
                    dataDictAnnot[rec.id]['End'] = dataDictAnnot[rec.id]['length']
                
                annotData = rec.letter_annotations["phred_quality"][int(dataDictAnnot[rec.id]['Begening']): int(dataDictAnnot[rec.id]['End'])]
                records[i].letter_annotations.pop("phred_quality")
                records[i].seq = records[i].seq[int(dataDictAnnot[rec.id]['Begening']): int(dataDictAnnot[rec.id]['End'])]
                records[i].letter_annotations["phred_quality"] = annotData
                newRecords.append(records[i])
            except KeyError:
                continue
        
        return newRecords
    
    
    for lines in fileData:
        if 'evalue' in lines:
            if ', Frame' in lines:
                eval = lines.split(' = ')[2].split(',')[0]
            else:
                eval = lines.split(' = ')[2]
            dataDict[lines.split(' : ')[0]] = (eval)
    
    recordsFasta = seqTrim(recordsFasta)
    
    for rec in recordsFasta:
        try:
            if float(1e-10) >= float(dataDict[rec.id]) > float(1e-30):
                recordsF1.append(rec)
            elif float(1e-30) >= float(dataDict[rec.id]) > float(1e-50):
                recordsF2.append(rec)
            if float(1e-50) >= float(dataDict[rec.id]):
                recordsF3.append(rec)
        
        except KeyError:
            continue
    
    for rec in recordsQual:
        if rec.id in dataDict.keys():
            if float(1e-10) >= float(dataDict[rec.id]) > float(1e-30):
                recordsQ1.append(rec)
            elif float(1e-30) >= float(dataDict[rec.id]) > float(1e-50):
                recordsQ2.append(rec)
            if float(1e-50) >= float(dataDict[rec.id]):
                recordsQ3.append(rec)
    
    recordsQ1 = qualEdit(recordsQ1)
    recordsQ2 = qualEdit(recordsQ2)
    recordsQ3 = qualEdit(recordsQ3)
    
    try:
        os.mkdir(outputFolder)
    except IOError:
        pass
    
    with open(outputFolder + '/seq_e_val_10-30.qual', 'w') as fp, open(outputFolder + '/seq_e_val_30-50.qual', 'w') as fq, open(outputFolder + '/seq_e_val_50-last.qual', 'w') as fr:
        SeqIO.write(recordsQ1, fp, 'qual')
        SeqIO.write(recordsQ2, fq, 'qual')
        SeqIO.write(recordsQ3, fr, 'qual')
    
    with open(outputFolder + '/seq_e_val_10-30.fasta', 'w') as fp, open(outputFolder + '/seq_e_val_30-50.fasta', 'w') as fq, open(outputFolder + '/seq_e_val_50-last.fasta', 'w') as fr:
        SeqIO.write(recordsF1, fp, 'fasta')
        SeqIO.write(recordsF2, fq, 'fasta')
        SeqIO.write(recordsF3, fr, 'fasta')
    
    
    files = ['Sequences.fasta', 'RawQual.qual', 'sequenceTagCheck.txt', 'Annotations.txt', 'Plot_pre.png', 'Plot_post.png']
    
    for file in files:
        shutil.copy2(files[i], outputFolder)
    
    print("All Done\n")



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
        handleQ = open('RawQual.qual', 'rU'); recordsQ = SeqIO.to_dict(SeqIO.parse(handleQ, 'qual'))
    
    elif args.fa != None:
        """Executes if Fasta input file supplied"""
        handle = open(args.fa, 'rU'); records = list(SeqIO.parse(handle, 'fasta'))
    
    handleF = open(args.frd, 'rU'); recordF = list(SeqIO.parse(handleF, 'fasta'))
    handleR = open(args.rev, 'rU'); recordR = list(SeqIO.parse(handleR, 'fasta'))
    
    if "/" not in records[0].id:
        for i, rec in enumerate(records):
            records[i].id = records[i].id + '/' + str(i)
        with open(args.fa.split('.')[0] + '_Edited.txt', 'w') as fp:
            SeqIO.write(records, fp, 'fasta')
        print("Your input sequence IDs has been changed. Please check %s_Edited.txt file for new sequence IDs\n\n" %args.fa.split('.')[0])
    
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
        print("Running analysis for -> %s\n"%val.id)
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
                records[i].seq = val.seq[endPosListR[0] + 1: -25]
            else:
                records[i].seq = val.seq[24: startPosListR[-1]]
        elif fFlag == False and rFlag == True:
            if endPosListF[0]/len(val.seq) < 0.7:
                records[i].seq = val.seq[endPosListF[0] + 1: -25]
            else:
                records[i].seq = val.seq[24: startPosListF[-1]]
        elif fFlag == True and rFlag == True:
            records[i].seq = val.seq[24: -25]
        else:
            pass
    
    print("Primer mapping done!!\n")
    
    with open('BlastInput.fas', 'w') as fp, open('FastqOutput.fq', 'w') as fq:
        SeqIO.write(records, fp, 'fasta')
        #SeqIO.write(recordsQ, fq, 'qual')
    
    with open('Annotations.txt', 'w') as fp:
        for val in records:
            fp.write('%s: %s\t\t%s\t\t%s\t\t%s\n' %(val.id, val.annotations, posDict[val.id][0], posDict[val.id][1], posDict[val.id][2]))
    
    
    print("Initiating second stage blastx search\n")
    blastxOR('BlastInput.fas')
    print("All Done. OR tags are stored in sequenceTag.txt file\n")
    
    
    evalPre = []; evalPost = []
    dataPre = open('sequenceTagCheck.txt', 'r').readlines(); dataPost = open('sequenceTag.txt', 'r').readlines()
    for linesPre, linesPost in zip(dataPre, dataPost):
        if 'evalue' in linesPre:
            evalPre.append(float(linesPre.split(',')[1].split(' = ')[1]))
        if 'evalue' in linesPost:
            evalPost.append(float(linesPost.split(',')[1].split(' = ')[1]))

    preDict = defaultdict(list); postDict = defaultdict(list)
    dictionaries = {'pre': preDict, 'post': postDict}

    preDict['10-20'] = 0; preDict['20-30'] = 0; preDict['30-40'] = 0; preDict['40-50'] = 0;
    preDict['50-60'] = 0; preDict['60-70'] = 0; preDict['70-80'] = 0; preDict['80-90'] = 0;
    preDict['90-100'] = 0; preDict['100-110'] = 0; preDict['110-120'] = 0; preDict['120-130'] = 0;
    preDict['130-140'] = 0; preDict['140-150'] = 0; preDict['150-160'] = 0; preDict['160-170'] = 0;
    preDict['170-180'] = 0; preDict['180-190'] = 0; preDict['190-200'] = 0;

    postDict['10-20'] = 0; postDict['20-30'] = 0; postDict['30-40'] = 0; postDict['40-50'] = 0;
    postDict['50-60'] = 0; postDict['60-70'] = 0; postDict['70-80'] = 0; postDict['80-90'] = 0;
    postDict['90-100'] = 0; postDict['100-110'] = 0; postDict['110-120'] = 0; postDict['120-130'] = 0;
    postDict['130-140'] = 0; postDict['140-150'] = 0; postDict['150-160'] = 0; postDict['160-170'] = 0;
    postDict['170-180'] = 0; postDict['180-190'] = 0; postDict['190-200'] = 0;

    with open('EvalPre.txt', 'w') as fp, open('EvalPost.txt', 'w') as fq:
        for val1, val2 in zip(evalPre, evalPost):
            fp.write('%s\n' %val1); fq.write('%s\n' %val2)
            preDict['10-20'] = preDict['10-20'] + 1 if val1 <= float(1e-10) and val1 > float(1e-20) else preDict['10-20']
            preDict['20-30'] = preDict['20-30'] + 1 if val1 <= float(1e-20) and val1 > float(1e-30) else preDict['20-30']
            preDict['30-40'] = preDict['30-40'] + 1 if val1 <= float(1e-30) and val1 > float(1e-40) else preDict['30-40']
            preDict['40-50'] = preDict['40-50'] + 1 if val1 <= float(1e-40) and val1 > float(1e-50) else preDict['40-50']
            preDict['50-60'] = preDict['50-60'] + 1 if val1 <= float(1e-50) and val1 > float(1e-60) else preDict['50-60']
            preDict['60-70'] = preDict['60-70'] + 1 if val1 <= float(1e-60) and val1 > float(1e-70) else preDict['60-70']
            preDict['70-80'] = preDict['70-80'] + 1 if val1 <= float(1e-70) and val1 > float(1e-80) else preDict['70-80']
            preDict['80-90'] = preDict['80-90'] + 1 if val1 <= float(1e-80) and val1 > float(1e-90) else preDict['80-90']
            preDict['90-100'] = preDict['90-100'] + 1 if val1 <= float(1e-90) and val1 > float(1e-100) else preDict['90-100']
            preDict['100-110'] = preDict['100-110'] + 1 if val1 <= float(1e-100) and val1 > float(1e-110) else preDict['100-110']
            preDict['110-120'] = preDict['110-120'] + 1 if val1 <= float(1e-110) and val1 > float(1e-120) else preDict['110-120']
            preDict['120-130'] = preDict['120-130'] + 1 if val1 <= float(1e-120) and val1 > float(1e-130) else preDict['120-130']
            preDict['130-140'] = preDict['130-140'] + 1 if val1 <= float(1e-130) and val1 > float(1e-140) else preDict['130-140']
            preDict['140-150'] = preDict['140-150'] + 1 if val1 <= float(1e-140) and val1 > float(1e-150) else preDict['140-150']
            preDict['150-160'] = preDict['150-160'] + 1 if val1 <= float(1e-150) and val1 > float(1e-160) else preDict['150-160']
            preDict['160-170'] = preDict['160-170'] + 1 if val1 <= float(1e-160) and val1 > float(1e-170) else preDict['160-170']
            preDict['160-170'] = preDict['160-170'] + 1 if val1 <= float(1e-160) and val1 > float(1e-170) else preDict['160-170']
            preDict['170-180'] = preDict['170-180'] + 1 if val1 <= float(1e-170) and val1 > float(1e-180) else preDict['170-180']
            preDict['180-190'] = preDict['180-190'] + 1 if val1 <= float(1e-180) and val1 > float(1e-190) else preDict['180-190']
            preDict['190-200'] = preDict['190-200'] + 1 if val1 <= float(1e-190) and val1 > float(1e-200) else preDict['190-200']
        
            postDict['10-20'] = postDict['10-20'] + 1 if val2 <= float(1e-10) and val2 > float(1e-20) else postDict['10-20']
            postDict['20-30'] = postDict['20-30'] + 1 if val2 <= float(1e-20) and val2 > float(1e-30) else postDict['20-30']
            postDict['30-40'] = postDict['30-40'] + 1 if val2 <= float(1e-30) and val2 > float(1e-40) else postDict['30-40']
            postDict['40-50'] = postDict['40-50'] + 1 if val2 <= float(1e-40) and val2 > float(1e-50) else postDict['40-50']
            postDict['50-60'] = postDict['50-60'] + 1 if val2 <= float(1e-50) and val2 > float(1e-60) else postDict['50-60']
            postDict['60-70'] = postDict['60-70'] + 1 if val2 <= float(1e-60) and val2 > float(1e-70) else postDict['60-70']
            postDict['70-80'] = postDict['70-80'] + 1 if val2 <= float(1e-70) and val2 > float(1e-80) else postDict['70-80']
            postDict['80-90'] = postDict['80-90'] + 1 if val2 <= float(1e-80) and val2 > float(1e-90) else postDict['80-90']
            postDict['90-100'] = postDict['90-100'] + 1 if val2 <= float(1e-90) and val2 > float(1e-100) else postDict['90-100']
            postDict['100-110'] = postDict['100-110'] + 1 if val2 <= float(1e-100) and val2 > float(1e-110) else postDict['100-110']
            postDict['110-120'] = postDict['110-120'] + 1 if val2 <= float(1e-110) and val2 > float(1e-120) else postDict['110-120']
            postDict['120-130'] = postDict['120-130'] + 1 if val2 <= float(1e-120) and val2 > float(1e-130) else postDict['120-130']
            postDict['130-140'] = postDict['130-140'] + 1 if val2 <= float(1e-130) and val2 > float(1e-140) else postDict['130-140']
            postDict['140-150'] = postDict['140-150'] + 1 if val2 <= float(1e-140) and val2 > float(1e-150) else postDict['140-150']
            postDict['150-160'] = postDict['150-160'] + 1 if val2 <= float(1e-150) and val2 > float(1e-160) else postDict['150-160']
            postDict['160-170'] = postDict['160-170'] + 1 if val2 <= float(1e-160) and val2 > float(1e-170) else postDict['160-170']
            postDict['160-170'] = postDict['160-170'] + 1 if val2 <= float(1e-160) and val2 > float(1e-170) else postDict['160-170']
            postDict['170-180'] = postDict['170-180'] + 1 if val2 <= float(1e-170) and val2 > float(1e-180) else postDict['170-180']
            postDict['180-190'] = postDict['180-190'] + 1 if val2 <= float(1e-180) and val2 > float(1e-190) else postDict['180-190']
            postDict['190-200'] = postDict['190-200'] + 1 if val2 <= float(1e-190) and val2 > float(1e-200) else postDict['190-200']


    for key, val in dictionaries.items():
        total = 0
        for inkey, inval in val.items():
            total = total + inval
        print key, total

    def dictUpdate(dictData):
        for key, val in dictData.items():
            dictData[int(key.split('-')[0])] = dictData.pop(key)
        dictData = collections.OrderedDict(sorted(dictData.items()))
        return dictData

    preDict = dictUpdate(preDict)
    postDict = dictUpdate(postDict)

    for key, d in dictionaries.items():
        fileName = "Plot_%s.png" %(key)
        d = collections.OrderedDict(sorted(d.items()))
        X = np.arange(len(d))
        pl.bar(X, d.values(), align='center', width=0.5)
        pl.xticks(X, sorted(d.keys()))
        pl.xticks(fontsize=8, rotation=90)
        ymax = max(d.values()) + 1
        pl.ylim(0, ymax)
        pl.savefig(fileName, format='png')
        pl.clf()

    processOut(args.o)


if __name__ == "__main__":
    main()







