#!/usr/bin/env python3

from __future__ import print_function
import pysam
import argparse

import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class NameEquivalenciesContainer:
    def __init__(self):
        self.map = {}

    def addEquivalency(self, name1, name2):
        self.map[name1] = name2
        self.map[name2] = name1

    def __contains__(self, item):
        return item in self.map

    def __getitem__(self, item):
        return self.map[item]


name_equivalencies = NameEquivalenciesContainer()
name_equivalencies.addEquivalency('chr1', '1')
name_equivalencies.addEquivalency('chr2', '2')
name_equivalencies.addEquivalency('chr3', '3')
name_equivalencies.addEquivalency('chr4', '4')
name_equivalencies.addEquivalency('chr5', '5')
name_equivalencies.addEquivalency('chr6', '6')
name_equivalencies.addEquivalency('chr7', '7')
name_equivalencies.addEquivalency('chr8', '8')
name_equivalencies.addEquivalency('chr9', '9')
name_equivalencies.addEquivalency('chr10', '10')
name_equivalencies.addEquivalency('chr11', '11')
name_equivalencies.addEquivalency('chr12', '12')
name_equivalencies.addEquivalency('chr13', '13')
name_equivalencies.addEquivalency('chr14', '14')
name_equivalencies.addEquivalency('chr15', '15')
name_equivalencies.addEquivalency('chr16', '16')
name_equivalencies.addEquivalency('chr17', '17')
name_equivalencies.addEquivalency('chr18', '18')
name_equivalencies.addEquivalency('chr19', '19')
name_equivalencies.addEquivalency('chr20', '20')
name_equivalencies.addEquivalency('chr21', '21')
name_equivalencies.addEquivalency('chr22', '22')
name_equivalencies.addEquivalency('chrX', 'X')
name_equivalencies.addEquivalency('chrY', 'Y')
name_equivalencies.addEquivalency('chrM', 'M')
name_equivalencies.addEquivalency('chr1_KI270706v1_random', '1_KI270706v1_random')
name_equivalencies.addEquivalency('chr1_KI270707v1_random', '1_KI270707v1_random')
name_equivalencies.addEquivalency('chr1_KI270708v1_random', '1_KI270708v1_random')
name_equivalencies.addEquivalency('chr1_KI270709v1_random', '1_KI270709v1_random')
name_equivalencies.addEquivalency('chr1_KI270710v1_random', '1_KI270710v1_random')
name_equivalencies.addEquivalency('chr1_KI270711v1_random', '1_KI270711v1_random')
name_equivalencies.addEquivalency('chr1_KI270712v1_random', '1_KI270712v1_random')
name_equivalencies.addEquivalency('chr1_KI270713v1_random', '1_KI270713v1_random')
name_equivalencies.addEquivalency('chr1_KI270714v1_random', '1_KI270714v1_random')
name_equivalencies.addEquivalency('chr2_KI270715v1_random', '2_KI270715v1_random')
name_equivalencies.addEquivalency('chr2_KI270716v1_random', '2_KI270716v1_random')
name_equivalencies.addEquivalency('chr3_GL000221v1_random', '3_GL000221v1_random')
name_equivalencies.addEquivalency('chr4_GL000008v2_random', '4_GL000008v2_random')
name_equivalencies.addEquivalency('chr5_GL000208v1_random', '5_GL000208v1_random')
name_equivalencies.addEquivalency('chr9_KI270717v1_random', '9_KI270717v1_random')
name_equivalencies.addEquivalency('chr9_KI270718v1_random', '9_KI270718v1_random')
name_equivalencies.addEquivalency('chr9_KI270719v1_random', '9_KI270719v1_random')
name_equivalencies.addEquivalency('chr9_KI270720v1_random', '9_KI270720v1_random')
name_equivalencies.addEquivalency('chr11_KI270721v1_random', '11_KI270721v1_random')
name_equivalencies.addEquivalency('chr14_GL000009v2_random', '14_GL000009v2_random')
name_equivalencies.addEquivalency('chr14_GL000225v1_random', '14_GL000225v1_random')
name_equivalencies.addEquivalency('chr14_KI270722v1_random', '14_KI270722v1_random')
name_equivalencies.addEquivalency('chr14_GL000194v1_random', '14_GL000194v1_random')
name_equivalencies.addEquivalency('chr14_KI270723v1_random', '14_KI270723v1_random')
name_equivalencies.addEquivalency('chr14_KI270724v1_random', '14_KI270724v1_random')
name_equivalencies.addEquivalency('chr14_KI270725v1_random', '14_KI270725v1_random')
name_equivalencies.addEquivalency('chr14_KI270726v1_random', '14_KI270726v1_random')
name_equivalencies.addEquivalency('chr15_KI270727v1_random', '15_KI270727v1_random')
name_equivalencies.addEquivalency('chr16_KI270728v1_random', '16_KI270728v1_random')
name_equivalencies.addEquivalency('chr17_GL000205v2_random', '17_GL000205v2_random')
name_equivalencies.addEquivalency('chr17_KI270729v1_random', '17_KI270729v1_random')
name_equivalencies.addEquivalency('chr17_KI270730v1_random', '17_KI270730v1_random')
name_equivalencies.addEquivalency('chr22_KI270731v1_random', '22_KI270731v1_random')
name_equivalencies.addEquivalency('chr22_KI270732v1_random', '22_KI270732v1_random')
name_equivalencies.addEquivalency('chr22_KI270733v1_random', '22_KI270733v1_random')
name_equivalencies.addEquivalency('chr22_KI270734v1_random', '22_KI270734v1_random')
name_equivalencies.addEquivalency('chr22_KI270735v1_random', '22_KI270735v1_random')
name_equivalencies.addEquivalency('chr22_KI270736v1_random', '22_KI270736v1_random')
name_equivalencies.addEquivalency('chr22_KI270737v1_random', '22_KI270737v1_random')
name_equivalencies.addEquivalency('chr22_KI270738v1_random', '22_KI270738v1_random')
name_equivalencies.addEquivalency('chr22_KI270739v1_random', '22_KI270739v1_random')
name_equivalencies.addEquivalency('chrY_KI270740v1_random', 'Y_KI270740v1_random')
name_equivalencies.addEquivalency('chrUn_KI270302v1', 'Un_KI270302v1')
name_equivalencies.addEquivalency('chrUn_KI270304v1', 'Un_KI270304v1')
name_equivalencies.addEquivalency('chrUn_KI270303v1', 'Un_KI270303v1')
name_equivalencies.addEquivalency('chrUn_KI270305v1', 'Un_KI270305v1')
name_equivalencies.addEquivalency('chrUn_KI270322v1', 'Un_KI270322v1')
name_equivalencies.addEquivalency('chrUn_KI270320v1', 'Un_KI270320v1')
name_equivalencies.addEquivalency('chrUn_KI270310v1', 'Un_KI270310v1')
name_equivalencies.addEquivalency('chrUn_KI270316v1', 'Un_KI270316v1')
name_equivalencies.addEquivalency('chrUn_KI270315v1', 'Un_KI270315v1')
name_equivalencies.addEquivalency('chrUn_KI270312v1', 'Un_KI270312v1')
name_equivalencies.addEquivalency('chrUn_KI270311v1', 'Un_KI270311v1')
name_equivalencies.addEquivalency('chrUn_KI270317v1', 'Un_KI270317v1')
name_equivalencies.addEquivalency('chrUn_KI270412v1', 'Un_KI270412v1')
name_equivalencies.addEquivalency('chrUn_KI270411v1', 'Un_KI270411v1')
name_equivalencies.addEquivalency('chrUn_KI270414v1', 'Un_KI270414v1')
name_equivalencies.addEquivalency('chrUn_KI270419v1', 'Un_KI270419v1')
name_equivalencies.addEquivalency('chrUn_KI270418v1', 'Un_KI270418v1')
name_equivalencies.addEquivalency('chrUn_KI270420v1', 'Un_KI270420v1')
name_equivalencies.addEquivalency('chrUn_KI270424v1', 'Un_KI270424v1')
name_equivalencies.addEquivalency('chrUn_KI270417v1', 'Un_KI270417v1')
name_equivalencies.addEquivalency('chrUn_KI270422v1', 'Un_KI270422v1')
name_equivalencies.addEquivalency('chrUn_KI270423v1', 'Un_KI270423v1')
name_equivalencies.addEquivalency('chrUn_KI270425v1', 'Un_KI270425v1')
name_equivalencies.addEquivalency('chrUn_KI270429v1', 'Un_KI270429v1')
name_equivalencies.addEquivalency('chrUn_KI270442v1', 'Un_KI270442v1')
name_equivalencies.addEquivalency('chrUn_KI270466v1', 'Un_KI270466v1')
name_equivalencies.addEquivalency('chrUn_KI270465v1', 'Un_KI270465v1')
name_equivalencies.addEquivalency('chrUn_KI270467v1', 'Un_KI270467v1')
name_equivalencies.addEquivalency('chrUn_KI270435v1', 'Un_KI270435v1')
name_equivalencies.addEquivalency('chrUn_KI270438v1', 'Un_KI270438v1')
name_equivalencies.addEquivalency('chrUn_KI270468v1', 'Un_KI270468v1')
name_equivalencies.addEquivalency('chrUn_KI270510v1', 'Un_KI270510v1')
name_equivalencies.addEquivalency('chrUn_KI270509v1', 'Un_KI270509v1')
name_equivalencies.addEquivalency('chrUn_KI270518v1', 'Un_KI270518v1')
name_equivalencies.addEquivalency('chrUn_KI270508v1', 'Un_KI270508v1')
name_equivalencies.addEquivalency('chrUn_KI270516v1', 'Un_KI270516v1')
name_equivalencies.addEquivalency('chrUn_KI270512v1', 'Un_KI270512v1')
name_equivalencies.addEquivalency('chrUn_KI270519v1', 'Un_KI270519v1')
name_equivalencies.addEquivalency('chrUn_KI270522v1', 'Un_KI270522v1')
name_equivalencies.addEquivalency('chrUn_KI270511v1', 'Un_KI270511v1')
name_equivalencies.addEquivalency('chrUn_KI270515v1', 'Un_KI270515v1')
name_equivalencies.addEquivalency('chrUn_KI270507v1', 'Un_KI270507v1')
name_equivalencies.addEquivalency('chrUn_KI270517v1', 'Un_KI270517v1')
name_equivalencies.addEquivalency('chrUn_KI270529v1', 'Un_KI270529v1')
name_equivalencies.addEquivalency('chrUn_KI270528v1', 'Un_KI270528v1')
name_equivalencies.addEquivalency('chrUn_KI270530v1', 'Un_KI270530v1')
name_equivalencies.addEquivalency('chrUn_KI270539v1', 'Un_KI270539v1')
name_equivalencies.addEquivalency('chrUn_KI270538v1', 'Un_KI270538v1')
name_equivalencies.addEquivalency('chrUn_KI270544v1', 'Un_KI270544v1')
name_equivalencies.addEquivalency('chrUn_KI270548v1', 'Un_KI270548v1')
name_equivalencies.addEquivalency('chrUn_KI270583v1', 'Un_KI270583v1')
name_equivalencies.addEquivalency('chrUn_KI270587v1', 'Un_KI270587v1')
name_equivalencies.addEquivalency('chrUn_KI270580v1', 'Un_KI270580v1')
name_equivalencies.addEquivalency('chrUn_KI270581v1', 'Un_KI270581v1')
name_equivalencies.addEquivalency('chrUn_KI270579v1', 'Un_KI270579v1')
name_equivalencies.addEquivalency('chrUn_KI270589v1', 'Un_KI270589v1')
name_equivalencies.addEquivalency('chrUn_KI270590v1', 'Un_KI270590v1')
name_equivalencies.addEquivalency('chrUn_KI270584v1', 'Un_KI270584v1')
name_equivalencies.addEquivalency('chrUn_KI270582v1', 'Un_KI270582v1')
name_equivalencies.addEquivalency('chrUn_KI270588v1', 'Un_KI270588v1')
name_equivalencies.addEquivalency('chrUn_KI270593v1', 'Un_KI270593v1')
name_equivalencies.addEquivalency('chrUn_KI270591v1', 'Un_KI270591v1')
name_equivalencies.addEquivalency('chrUn_KI270330v1', 'Un_KI270330v1')
name_equivalencies.addEquivalency('chrUn_KI270329v1', 'Un_KI270329v1')
name_equivalencies.addEquivalency('chrUn_KI270334v1', 'Un_KI270334v1')
name_equivalencies.addEquivalency('chrUn_KI270333v1', 'Un_KI270333v1')
name_equivalencies.addEquivalency('chrUn_KI270335v1', 'Un_KI270335v1')
name_equivalencies.addEquivalency('chrUn_KI270338v1', 'Un_KI270338v1')
name_equivalencies.addEquivalency('chrUn_KI270340v1', 'Un_KI270340v1')
name_equivalencies.addEquivalency('chrUn_KI270336v1', 'Un_KI270336v1')
name_equivalencies.addEquivalency('chrUn_KI270337v1', 'Un_KI270337v1')
name_equivalencies.addEquivalency('chrUn_KI270363v1', 'Un_KI270363v1')
name_equivalencies.addEquivalency('chrUn_KI270364v1', 'Un_KI270364v1')
name_equivalencies.addEquivalency('chrUn_KI270362v1', 'Un_KI270362v1')
name_equivalencies.addEquivalency('chrUn_KI270366v1', 'Un_KI270366v1')
name_equivalencies.addEquivalency('chrUn_KI270378v1', 'Un_KI270378v1')
name_equivalencies.addEquivalency('chrUn_KI270379v1', 'Un_KI270379v1')
name_equivalencies.addEquivalency('chrUn_KI270389v1', 'Un_KI270389v1')
name_equivalencies.addEquivalency('chrUn_KI270390v1', 'Un_KI270390v1')
name_equivalencies.addEquivalency('chrUn_KI270387v1', 'Un_KI270387v1')
name_equivalencies.addEquivalency('chrUn_KI270395v1', 'Un_KI270395v1')
name_equivalencies.addEquivalency('chrUn_KI270396v1', 'Un_KI270396v1')
name_equivalencies.addEquivalency('chrUn_KI270388v1', 'Un_KI270388v1')
name_equivalencies.addEquivalency('chrUn_KI270394v1', 'Un_KI270394v1')
name_equivalencies.addEquivalency('chrUn_KI270386v1', 'Un_KI270386v1')
name_equivalencies.addEquivalency('chrUn_KI270391v1', 'Un_KI270391v1')
name_equivalencies.addEquivalency('chrUn_KI270383v1', 'Un_KI270383v1')
name_equivalencies.addEquivalency('chrUn_KI270393v1', 'Un_KI270393v1')
name_equivalencies.addEquivalency('chrUn_KI270384v1', 'Un_KI270384v1')
name_equivalencies.addEquivalency('chrUn_KI270392v1', 'Un_KI270392v1')
name_equivalencies.addEquivalency('chrUn_KI270381v1', 'Un_KI270381v1')
name_equivalencies.addEquivalency('chrUn_KI270385v1', 'Un_KI270385v1')
name_equivalencies.addEquivalency('chrUn_KI270382v1', 'Un_KI270382v1')
name_equivalencies.addEquivalency('chrUn_KI270376v1', 'Un_KI270376v1')
name_equivalencies.addEquivalency('chrUn_KI270374v1', 'Un_KI270374v1')
name_equivalencies.addEquivalency('chrUn_KI270372v1', 'Un_KI270372v1')
name_equivalencies.addEquivalency('chrUn_KI270373v1', 'Un_KI270373v1')
name_equivalencies.addEquivalency('chrUn_KI270375v1', 'Un_KI270375v1')
name_equivalencies.addEquivalency('chrUn_KI270371v1', 'Un_KI270371v1')
name_equivalencies.addEquivalency('chrUn_KI270448v1', 'Un_KI270448v1')
name_equivalencies.addEquivalency('chrUn_KI270521v1', 'Un_KI270521v1')
name_equivalencies.addEquivalency('chrUn_GL000195v1', 'Un_GL000195v1')
name_equivalencies.addEquivalency('chrUn_GL000219v1', 'Un_GL000219v1')
name_equivalencies.addEquivalency('chrUn_GL000220v1', 'Un_GL000220v1')
name_equivalencies.addEquivalency('chrUn_GL000224v1', 'Un_GL000224v1')
name_equivalencies.addEquivalency('chrUn_KI270741v1', 'Un_KI270741v1')
name_equivalencies.addEquivalency('chrUn_GL000226v1', 'Un_GL000226v1')
name_equivalencies.addEquivalency('chrUn_GL000213v1', 'Un_GL000213v1')
name_equivalencies.addEquivalency('chrUn_KI270743v1', 'Un_KI270743v1')
name_equivalencies.addEquivalency('chrUn_KI270744v1', 'Un_KI270744v1')
name_equivalencies.addEquivalency('chrUn_KI270745v1', 'Un_KI270745v1')
name_equivalencies.addEquivalency('chrUn_KI270746v1', 'Un_KI270746v1')
name_equivalencies.addEquivalency('chrUn_KI270747v1', 'Un_KI270747v1')
name_equivalencies.addEquivalency('chrUn_KI270748v1', 'Un_KI270748v1')
name_equivalencies.addEquivalency('chrUn_KI270749v1', 'Un_KI270749v1')
name_equivalencies.addEquivalency('chrUn_KI270750v1', 'Un_KI270750v1')
name_equivalencies.addEquivalency('chrUn_KI270751v1', 'Un_KI270751v1')
name_equivalencies.addEquivalency('chrUn_KI270752v1', 'Un_KI270752v1')
name_equivalencies.addEquivalency('chrUn_KI270753v1', 'Un_KI270753v1')
name_equivalencies.addEquivalency('chrUn_KI270754v1', 'Un_KI270754v1')
name_equivalencies.addEquivalency('chrUn_KI270755v1', 'Un_KI270755v1')
name_equivalencies.addEquivalency('chrUn_KI270756v1', 'Un_KI270756v1')
name_equivalencies.addEquivalency('chrUn_KI270757v1', 'Un_KI270757v1')
name_equivalencies.addEquivalency('chrUn_GL000214v1', 'Un_GL000214v1')
name_equivalencies.addEquivalency('chrUn_KI270742v1', 'Un_KI270742v1')
name_equivalencies.addEquivalency('chrUn_GL000216v2', 'Un_GL000216v2')
name_equivalencies.addEquivalency('chrUn_GL000218v1', 'Un_GL000218v1')
name_equivalencies.addEquivalency('chrEBV', 'EBV')

requiredSequences = [
    'chr1',
    'chr2',
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    'chr8',
    'chr9',
    'chr10',
    'chr11',
    'chr12',
    'chr13',
    'chr14',
    'chr15',
    'chr16',
    'chr17',
    'chr18',
    'chr19',
    'chr20',
    'chr21',
    'chr22',
    'chrX',
    'chrY'
]


def getSequenceDictFromBamFile(inputPath) -> dict:
    samfile = pysam.AlignmentFile(inputPath, "rb")
    header_dict = samfile.header.as_dict()
    result = {}
    sequences = header_dict['SQ']
    for sequence in sequences:
        sequenceName = sequence['SN']
        sequenceLengh = sequence['LN']
        if 'M5' in sequence:
            sequenceM5 = sequence['M5']
        else:
            sequenceM5 = None
        result[sequenceName] = (sequenceLengh, sequenceM5)
    return result


def getSequenceDictFromReferenceDict(inputPath) -> dict:
    with open(inputPath, "r") as referenceDict:
        result = {}
        for line in referenceDict:
            if line.startswith('@HD'):
                continue
            line = line[4:]
            elements = [element.strip() for element in line.split("\t")]

            sequenceName = None
            sequenceLength = None
            sequenceM5 = None

            for element in elements:
                if element.startswith("SN:"):
                    if sequenceName is not None:
                        exit(-1)
                    sequenceName = element[3:]
                elif element.startswith("LN:"):
                    if sequenceLength is not None:
                        exit(-1)
                    sequenceLength = int(element[3:])
                elif element.startswith("M5:"):
                    if sequenceM5 is not None:
                        exit(-1)
                    sequenceM5 = element[3:]
            result[sequenceName] = (sequenceLength, sequenceM5)
        return result


def compareSequencesDict(dict1: dict, dict2: dict, requiredSequences):
    keys_in_both = set()
    keys_missing_in_one = set()

    for sequenceName in dict1.keys():
        if sequenceName in dict2 or (sequenceName in name_equivalencies and name_equivalencies[sequenceName] in dict2):
            keys_in_both.add(sequenceName)
        else:
            keys_missing_in_one.add(sequenceName)

    for sequenceName in dict2.keys():
        if sequenceName in dict1 or (sequenceName in name_equivalencies and name_equivalencies[sequenceName] in dict1):
            keys_in_both.add(sequenceName)
        else:
            keys_missing_in_one.add(sequenceName)

    for sequenceName in keys_missing_in_one:
        if sequenceName in requiredSequences:
            eprint("missing required sequence: ",sequenceName)
            return 1

    for sequenceName in keys_in_both:
        valDict1 = dict1[sequenceName] if sequenceName in dict1 else dict1[name_equivalencies[sequenceName]]
        valDict2 = dict2[sequenceName] if sequenceName in dict2 else dict2[name_equivalencies[sequenceName]]
        if valDict1[1] is None or valDict2[1] is None:
            valDict1 = (valDict1[0], None)
            valDict2 = (valDict2[0], None)

        if valDict1 != valDict2 and sequenceName in requiredSequences:
            eprint("incorrect required sequence: ",sequenceName)
            return 1
    return 0


def main():
    ap = argparse.ArgumentParser(description="idenitfy reference genome used for aligning")

    ap.add_argument("-b", "--bamFile", required=True)
    ap.add_argument("-d", "--dictFile", required=True)

    args = vars(ap.parse_args())
    # Create a dictionary of the shell arguments
    dictBam = getSequenceDictFromBamFile(args['bamFile'])
    dictReferenceDict = getSequenceDictFromReferenceDict(args['dictFile'])
    

    eprint("starting to print comparison")
    print(compareSequencesDict(dictBam, dictReferenceDict, requiredSequences))



if __name__ == '__main__':
    main()

