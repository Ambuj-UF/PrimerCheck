################################################################################################################
# Extract sequence and qual data and store it in different bin range                                           #
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


from Bio import SeqIO

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

fileData = open('sequenceTag.txt', 'r').readlines()

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

with open('seq_e_val_10-30.qual', 'w') as fp, open('seq_e_val_30-50.qual', 'w') as fq, open('seq_e_val_50-last.qual', 'w') as fr:
    SeqIO.write(recordsQ1, fp, 'qual')
    SeqIO.write(recordsQ2, fq, 'qual')
    SeqIO.write(recordsQ3, fr, 'qual')

with open('seq_e_val_10-30.fasta', 'w') as fp, open('seq_e_val_30-50.fasta', 'w') as fq, open('seq_e_val_50-last.fasta', 'w') as fr:
    SeqIO.write(recordsF1, fp, 'fasta')
    SeqIO.write(recordsF2, fq, 'fasta')
    SeqIO.write(recordsF3, fr, 'fasta')

