from os import path, makedirs


def Revcomp(sequence):
    comp = ""
    for letter in sequence:
        comp += pairBasePair(letter)
    return comp[::-1]

def pairBasePair(bp):
    switch = {
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G'
    }
    return switch.get(bp)


seq = dict()
with open(snakemake.input[0],'r') as inFile:
    dataList = inFile.readlines()
if not path.exists(path.dirname(snakemake.output[0])):
    makedirs(path.dirname(snakemake.output[0]))
with open(snakemake.output[0], 'w') as outFile:
    for i in range(0,len(dataList),2):
        seq ={
               "ID": dataList[i].strip('\n'),
               "barcode" : dataList[i+1].strip('\n')
              }
        print(f"{seq['ID']}", file=outFile)
        print(f"{Revcomp(seq['barcode'][:6])}{Revcomp(seq['barcode'][6:])}", file=outFile)
