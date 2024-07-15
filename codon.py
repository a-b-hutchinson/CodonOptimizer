import ast
import json
import re
from time import monotonic, sleep

import Bio
import optipyzer
import pandas as pd
import requests
from Bio import Entrez, SeqIO
from Bio.Restriction.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from IPython.display import clear_output
from orffinder import orffinder


def checkPeptides(entry) -> dict:
    """
    Searches 'features' section of uniprotKB entry for a transit or signal 
    peptide feature, returns dict with which peptides are present and where
    they are located.

    Parameters
    ----------
    entry: dict
        dict from getUniprotInfo method. will contain features key.
    """
    
    transit_present = False
    signal_present = False
    removal_regions = []

    for i in entry['features']:
        if 'Transit' in i['type']:
            transit_present = True
        if 'Signal' in i['type']:
            signal_present = True
        removal_regions.append([i['location']['start']['value'], i['location']['end']['value']])

    results = {
        'transit_peptide': transit_present,
        'signal_peptide': signal_present,
        'removed_region': removal_regions
    }

    return results

def getUniprotInfo(id) -> dict:
    """
    Querries UniprotKB using entry ID. Returns dict of response items.

    Parameters
    ----------
    id: str
        UniprotKB entry ID e.g. P05084
    """

    params = {
        'query': id,
        'fields': ['xref_refseq', 'sequence', 'ft_signal', 'ft_transit']
    }

    print(f'Processing [{id}]')

    res = requests.get(f'https://rest.uniprot.org/uniprotkb/{id}.fasta')
    fasta = res.content.decode('utf-8')
    if fasta == "Error messages\nThe 'accession' value has invalid format. It should be a valid UniProtKB accession":
        print('Issue with uniprot ID. Try removing isoform.')
    uniprotValue = ''.join(fasta.split('\n')[1:])
    les_nam = fasta[re.search('\\|.*\\|', fasta).end():]
    gene_symbol = re.search('(GN=.*? )',fasta)[0][3:].strip()
    name = re.search(' [A-Za-z0-9s][^ ]*', les_nam)[0][1:]
    response = requests.get('https://rest.uniprot.org/uniprotkb/search', params=params)

    try:
        entry = ast.literal_eval(str(response.content.decode('utf-8')))['results'][0]
    except ValueError as e:
        print(e, '\nAttempting to reprocess...')
        entry = ast.literal_eval(str(json.loads(response.content.decode('utf-8'))))['results'][0]

    potential_ids = pd.json_normalize(entry)['uniProtKBCrossReferences'][0]
    sequence_length = entry['sequence']['length']

    if len(entry['features']) > 0:
        peptides = checkPeptides(entry)
        transit_peptide = peptides['transit_peptide']
        signal_peptide = peptides['signal_peptide']
        removed_regions = peptides['removed_region']
    else:
        transit_peptide, signal_peptide, removed_regions = False, False, pd.NA

    results = {
        'uniprotValue': uniprotValue,
        'potential_ids': potential_ids,
        'sequence_length': sequence_length,
        'gene_symbol': gene_symbol,
        'name': name,
        'transit_peptide': transit_peptide,
        'signal_peptide': signal_peptide,
        'removed_regions': removed_regions
    }
    return results

def matchAmino(uniprotValue, potential_ids) -> list:
    """
    Uses UniprotKB amino sequence and possible NCBI Accession IDs to find
    matching entry returns the matching NCBI ID, DNA sequence, uniprot value

    Parameters
    ----------
    uniprotValue: str
        uniprotKB entry ID
    
    potential_ids: dict
        potential_ids are provided from getUniprotInfo results
    """
    if potential_ids == None:
        print('No potential ids provided.')
        return 0
    
    pt_ids = [i['properties'][0]['value'] for i in potential_ids]
    potential_ids = [i.split('.')[0] for i in pt_ids] + potential_ids

    for count, test in enumerate(potential_ids):
        print(f'\nChecking {count+1} of {len(potential_ids)} potential matches')
        print(f'Trying to match [{test}]...')

        #ENTER EMAIL HERE
        Entrez.email = ''
        handle = Entrez.efetch(db='nuccore', id=test, rettype='gb', retmode='text')
        sequence = SeqIO.read(handle, "genbank")

        for feature in sequence.features:
            if feature.type=='CDS':
                ncbi_match = feature.qualifiers['translation'][0]

        if ncbi_match == uniprotValue:
            print(f'Match found with [{test}]')
            return [test, sequence, ncbi_match]
        elif count+1 == len(potential_ids):
            continue
        else:
            print('Not a match. Idling for 10 seconds...', end=' ')
            start = monotonic()
            while monotonic()-start < 10:
                print('\u2588', end=' ')
                sleep(1.0 - ((monotonic()-start) % 1.0))

def matchTranslation(orfs, uniprotValue, sequence) -> dict:
    """
    Takes dict of potential orf matches and uniprot Amino sequence. Compares
    cut orf to find matching sequence.

    Parameters
    ----------
    orfs : dict
        Dict of orfs site from orffinder utility

    uniprotValue : str
        Uniprot entry ID
    """

    matchedTranslation = {}
    
    for i in range(len(orfs)):
        if uniprotValue+'*' == str(sequence.seq[orfs[i]['start']-1:orfs[i]['end']-1].translate()):
            print(f'Matched ORF at {orfs[i]['start']-1} to {orfs[i]['end']}')
            matchedTranslation['sliced'] = str(sequence.seq[orfs[i]['start']-1:orfs[i]['end']-1])
            matchedTranslation['region'] = [orfs[i]['start']-1, orfs[i]['end']]

            return matchedTranslation
        
def optimizeSlice(slice) -> str:
    """
    Optimizing sequence for E. Coli.

    Parameters
    ----------
    slice: str
        get string of appropriate ORF section from matchTranslation
    """

    api = optipyzer.API()
    gblock = str(slice)

    result = api.optimize(
        seq=gblock,
        seq_type="dna",
        weights={"e_coli": 1}
    )

    optimized = result['optimized_sd']
    return optimized

def findSites(optimizedCodon) -> dict:
    """
    Parameters
    ----------
    optimizedCodon : str
        Provides E. Coli optimized codon from optimizedSlice(). Will return dictionary
        with whether enzyme sites are present or not.
    """

    results = {'ins':[], 'outs':[]}
    enzymes = ['BsaI', 'BbsI', 'BsmBI', 'Esp3I', 'SapI']
    batch = RestrictionBatch()
    for i in enzymes:
        batch.add(i)

    for i in enzymes:
        enzyme = batch.get(i)
        if enzyme.site in optimizedCodon:
            results['ins'].append(i)
            print(f'{i} in sequence')
        else:
            results['outs'].append(i)
            print(f'{i} NOT in sequence')

    return results

def removePeptides(sequence, removal_regions) -> str:
    """
    Parameters
    ----------
    sequence: str
    
    removal_regions: list
    """
    sequence = list(sequence)
    to_remove = []
    for i in removal_regions:
        try:
            to_remove = to_remove + [j for j in range(i[0], i[1]*3)]
        except TypeError as e:
            print(e,'\n', f'There was a problem trying to identify peptide start ({i[0]}) or stop ({i[1]})\nRebuilding sequence...')
            sequence = ''.join(sequence)
            return sequence
            
        to_remove = list(set(to_remove))
        to_remove.sort(reverse=True)

    for i in to_remove:
        sequence.pop(i-1)
    
    sequence = ''.join(sequence)
    sequence = 'ATG' + sequence

    return sequence

def main(df) -> pd.DataFrame:
    """
    Parameters
    ----------
    df : pd.DataFrame
        csv with uniprot_id column
    """
    
    uniprotId = []
    name = []
    transit_peptide = []
    signal_peptide = []
    removed_regions = []
    symbols = []
    uniprotValue = []
    sequence_length = []
    nih_id = []
    full_sequence = []
    cut_sequence = []
    amino_sequence = []
    orf = []
    optimized_codon = []
    enzymes = []
    

    for count, i in enumerate(df.uniprot_id):
        clear_output(wait=True)
        print(f'Checking {count+1} of {len(df.uniprot_id)}')
        uniprotInfo = getUniprotInfo(i)
        
        try:
            matchedId, sequence, amino = matchAmino(uniprotInfo['uniprotValue'], uniprotInfo['potential_ids'])
        except:
            uniprotValue.append(uniprotInfo['uniprotValue'])
            sequence_length.append(uniprotInfo['sequence_length'])
            nih_id.append(pd.NA)
            full_sequence.append(pd.NA)
            cut_sequence.append(pd.NA)
            amino_sequence.append(pd.NA)
            orf.append(pd.NA)
            optimized_codon.append(pd.NA)
            enzymes.append(pd.NA)
            continue
        orfs = orffinder.getORFs(sequence, remove_nested=True)
        res = matchTranslation(orfs, uniprotInfo['uniprotValue'], sequence)
        if uniprotInfo['transit_peptide'] or uniprotInfo['signal_peptide']:
            cut_seq = removePeptides(res['sliced'], uniprotInfo['removed_regions'])
            print(cut_seq)
        else:
            cut_seq = pd.NA
        try:
            codon = optimizeSlice(res['sliced'])
            enzyme_results = findSites(optimizedCodon=codon)
        except TypeError as e:
            res = {}
            res['sliced'] = pd.NA
            print('Something\'s wrong with the codon')
            print(e)
            sleep(2)
            codon = e
            enzyme_results = pd.NA
        
        uniprotId.append(i)
        name.append(uniprotInfo['name'])
        transit_peptide.append(uniprotInfo['transit_peptide'])
        signal_peptide.append(uniprotInfo['signal_peptide'])
        removed_regions.append(uniprotInfo['removed_regions'])
        symbols.append(uniprotInfo['gene_symbol'])
        uniprotValue.append(uniprotInfo['uniprotValue'])
        sequence_length.append(uniprotInfo['sequence_length'])
        nih_id.append(matchedId)
        full_sequence.append(str(sequence.seq))
        cut_sequence.append(cut_seq)
        amino_sequence.append(amino)
        orf.append(res['sliced'])
        optimized_codon.append(codon)
        enzymes.append(enzyme_results)


        

    result = pd.DataFrame([
        uniprotId,
        name,
        transit_peptide,
        signal_peptide,
        removed_regions,
        symbols,
        uniprotValue,
        sequence_length,
        nih_id,
        full_sequence,
        cut_sequence,
        amino_sequence,
        orf,
        optimized_codon,
        enzymes
    ], index=[
        'uniprotID',
        'name',
        'tr_peptide',
        'si_peptide',
        'peptide_regions',
        'gene',
        'up_sequence',
        'seq_length',
        'nih_id',
        'full_sequence',
        'cut_sequence',
        'amino_sequence',
        'orf',
        'optimized_codon',
        'enzymes'
    ])
    
    result=result.transpose()

    for i in result.index:
        if pd.isna(result.loc[i,'enzymes']):
            continue
        
        cols = result.loc[i, 'enzymes']['ins'] + result.loc[i, 'enzymes']['outs']
        cols.sort()

    for i in cols:
        result[i] = pd.NA

    for i in result.index:
        if pd.isna(result.loc[i, 'enzymes']):
            continue

        for j in result.loc[i, 'enzymes']['ins']:
            result.loc[i, j] = True

        for k in result.loc[i, 'enzymes']['outs']:
            result.loc[i, k] = False

    result = result.drop('enzymes', axis=1)
    
    clear_output(wait=True)
    print('Done')
    return result

df = pd.read_csv(input('Enter path to file containing uniprot IDs\n>'))
result = main(df)
result.to_csv('results.csv')