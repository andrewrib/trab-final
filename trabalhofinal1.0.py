## Importando bibliotecas 
import sys
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
entry = sys.argv[1::]
seq = str(entry[0].upper())
arqdes = str(entry[1])
banco = str(entry[2])
content = ''
cab = []
saida = r'C:\Users\Ok\Desktop\trabalho final progII\Saida.txt'
blastpath = r'C:\Program Files\NCBI\blast-2.11.0+\bin\blastx.exe'
file = open(arqdes)
hurdur = []
for linha in file:
    if '>' in linha:
        cab.append(linha)
        hurdur.append(content)
        content = ''
    else:
        content += linha
hurdur.append(content)
for i in range(len(hurdur)):
    if seq in hurdur[i]:
        aux = cab[i-1].split('\n')
        aux2 = aux[0] + '\n' + hurdur[i]
        nome = aux[0] + '.fasta'
        arq = open(nome[1::], 'w')
        arq.write(aux2)
        arq.close()
        blase = NcbiblastxCommandline(cmd=blastpath,
                                      query=nome[1::],
                                      subject=banco,
                                      evalue=0.05,
                                      outfmt=6,
                                      out=saida)
        stdout, stdeer = blase()
        result = pd.read_csv(saida,
                             sep='\t',
                             names=["Seq_ID",
                                    "Amno_ID",
                                    "pident",
                                    "length",
                                    "mismatch",
                                    "gapopen",
                                    "qstart",
                                    "qend",
                                    "sstart",
                                    "send",
                                    "E-value",
                                    "Bitscore"])
        score = result.sort_values('Bitscore', ascending=False)
        print('-' * 100)
        print(score.iloc[[0], [0, 1, 10, 11]])
file.close
