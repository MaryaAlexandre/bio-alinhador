from Bio import Entrez

Entrez.email = "maryaeduardaa333@gmail.com"


SEQUENCIAS = {
    "humano_insulina.fasta": "NM_001185097",     # Homo sapiens insulin (INS)
    "camundongo_insulina.fasta": "NM_008387",    # Mus musculus insulin II (Ins2)
}

def baixar_sequencia(accession, nome_arquivo):
    dados = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    conteudo = dados.read()
    dados.close()

    with open(nome_arquivo, "w") as f:
        f.write(conteudo)

    titulo = conteudo.split("\n")[0]
    print(titulo)
    print(f"Salvo em {nome_arquivo}\n")

for nome_arquivo, accession in SEQUENCIAS.items():
    baixar_sequencia(accession, nome_arquivo)
