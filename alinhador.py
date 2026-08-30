from Bio import SeqIO


def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n, m = len(seq1), len(seq2)

    matriz = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        matriz[i][0] = i * gap
    for j in range(m + 1):
        matriz[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = matriz[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            cima = matriz[i-1][j] + gap
            esquerda = matriz[i][j-1] + gap
            matriz[i][j] = max(diagonal, cima, esquerda)

    alinhado1, alinhado2 = "", ""
    i, j = n, m
    while i > 0 and j > 0:
        atual = matriz[i][j]
        diagonal = matriz[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)

        if atual == diagonal:
            alinhado1 = seq1[i-1] + alinhado1
            alinhado2 = seq2[j-1] + alinhado2
            i -= 1
            j -= 1
        elif atual == matriz[i-1][j] + gap:
            alinhado1 = seq1[i-1] + alinhado1
            alinhado2 = "-" + alinhado2
            i -= 1
        else:
            alinhado1 = "-" + alinhado1
            alinhado2 = seq2[j-1] + alinhado2
            j -= 1

    while i > 0:
        alinhado1 = seq1[i-1] + alinhado1
        alinhado2 = "-" + alinhado2
        i -= 1
    while j > 0:
        alinhado1 = "-" + alinhado1
        alinhado2 = seq2[j-1] + alinhado2
        j -= 1

    return alinhado1, alinhado2, matriz[n][m]


def gerar_html_alinhamento(alinhado1, alinhado2, largura_bloco=60):
    blocos_html = []
    for inicio in range(0, len(alinhado1), largura_bloco):
        bloco1 = alinhado1[inicio:inicio + largura_bloco]
        bloco2 = alinhado2[inicio:inicio + largura_bloco]

        html1, html2 = "", ""
        for c1, c2 in zip(bloco1, bloco2):
            if c1 == "-" or c2 == "-":
                cor = "#555555"
            elif c1 == c2:
                cor = "#2e7d32"
            else:
                cor = "#b71c1c"

            html1 += f'<span style="background-color:{cor}">{c1}</span>'
            html2 += f'<span style="background-color:{cor}">{c2}</span>'

        blocos_html.append(f'<div>{html1}</div><div>{html2}</div><div style="height:14px"></div>')

    corpo = "\n".join(blocos_html)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
body {{ font-family: monospace; font-size: 15px; background: #1e1e1e; color: #eee; padding: 24px; }}
</style>
</head>
<body>
{corpo}
</body>
</html>"""


registro_humano = SeqIO.read("humano_insulina.fasta", "fasta")
registro_camundongo = SeqIO.read("camundongo_insulina.fasta", "fasta")

seq_humano = str(registro_humano.seq)
seq_camundongo = str(registro_camundongo.seq)

print(f"Humano: {len(seq_humano)} bases")
print(f"Camundongo: {len(seq_camundongo)} bases")

alinhado1, alinhado2, score = needleman_wunsch(seq_humano, seq_camundongo)

print("Pontuação final:", score)

matches = sum(1 for a, b in zip(alinhado1, alinhado2) if a == b and a != "-")
identidade = matches / len(alinhado1) * 100
print(f"Identidade: {identidade:.1f}%")

with open("resultado_alinhamento.txt", "w") as f:
    f.write(alinhado1 + "\n" + alinhado2 + "\n")

print("Alinhamento salvo em resultado_alinhamento.txt")

html = gerar_html_alinhamento(alinhado1, alinhado2)
with open("alinhamento.html", "w") as f:
    f.write(html)
print("Visualização salva em alinhamento.html")