## O programa deve ser executado no terminal, de preferencia no miniconda.
## Utilize o 'python trabalhofinal2.py -h' para ler a descrição.

## Importando bibliotecas
import sys
import pandas as pd
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline

parser = argparse.ArgumentParser(description = 'Programa que realiza busca de um motivo de DNA de sequências desconhecidas \n OBS: É necessário trocar o caminho da saída e trocar o caminho do blastpath')
parser.add_argument('seq', type=str, help= 'Sequência escolhida pelo usuário para realizar a busca')
parser.add_argument('arqdes', type=str, help= 'Sequência desconhecida do usuário')
parser.add_argument('banco', type=str, help= 'Arquivo onde será comparado as sequências')
args = parser.parse_args()

content = '' ## criando variável;
cab = [] ## criando uma lista;
saida = r'C:\Users\Ok\Desktop\trabalho final progII\Saida.txt'
blastpath = r'C:\Program Files\NCBI\blast-2.11.0+\bin\blastx.exe'
file = open(args.arqdes) ## abrindo o arquivo desconhecido;
hurdur = [] ## criando uma lista;
for linha in file:
    if '>' in linha:
        cab.append(linha) ## pegando o nome do gene;
        hurdur.append(content) ## pega a sequência de DNA do arquivo e jogando na lista;
        content = '' ##limpando a variável;
    else:
        content += linha
hurdur.append(content) ## pegando o última sequência;
for i in range(len(hurdur)):
    if args.seq in hurdur[i]: ## realizando a busca;
        aux = cab[i-1].split('\n') ## separa os arquivos de acordo com a busca;
        aux2 = aux[0] + '\n' + hurdur[i] ## pegando o conteúdo;
        nome = aux[0] + '.fasta' ## criando o nome do arquivo;
        arq = open(nome[1::], 'w') ## criando o arquivo;
        arq.write(aux2) ## escrevendo no arquivo;
        arq.close() ## fechando o arquivo;
        blase = NcbiblastxCommandline(cmd=blastpath,
                                      query=nome[1::], ## o que quero buscar;
                                      subject=args.banco, ## Onde quero buscar;
                                      evalue=0.05, ## e-valor;
                                      outfmt=6, ## formato do arquivo de saida;
                                      out=saida) ## caminho da saida;
        stdout, stdeer = blase()
        result = pd.read_csv(saida,
                             sep='\t',
                             names=["Seq_ID", ## query (e.g., unknown gene) sequence id;
                                    "Amno_ID", ## subject (e.g., reference genome) sequence id;
                                    "pident", ## percentage of identical matches;
                                    "length", ## alignment length (sequence overlap);
                                    "mismatch", ## number of mismatches;
                                    "gapopen", ## number of gap openings;
                                    "qstart", ## start of alignment in query;
                                    "qend", ## end of alignment in query;
                                    "sstart", ## start of alignment in subject;
                                    "send", ## end of alignment in subject;
                                    "E-value", ## expect value;
                                    "Bitscore"]) ## bit score;
        score = result.sort_values('Bitscore', ascending=False)
        print('-' * 100) ## divisor;
        print(score.iloc[[0], [0, 1, 10, 11]]) ## printando os dados;
file.close ## fechando o programa.
