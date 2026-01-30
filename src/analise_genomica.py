from Bio import SeqIO
from collections import Counter
import os


def analisar_genoma_fago(genbank_file, output_dir):
    record = SeqIO.read(genbank_file, "genbank")

    phage_name = record.annotations.get("organism", "Unknown phage")
    accession = record.id
    genome_size = len(record.seq)

    cds_features = [f for f in record.features if f.type == "CDS"]
    total_cds = len(cds_features)

    categorias = {
        "estrutural": ["capsid", "head", "tail", "portal", "baseplate"],
        "replicacao": ["polymerase", "helicase", "primase", "ligase"],
        "lise": ["holin", "endolysin", "lysin"],
        "regulacao": ["repressor", "regulatory", "transcription"],
        "hipotetico": ["hypothetical"]
    }

    contagem_funcional = Counter()
    genes_info = []

    for cds in cds_features:
        produto = cds.qualifiers.get("product", ["unknown"])[0].lower()
        gene = cds.qualifiers.get("gene", ["-"])[0]
        locus = cds.qualifiers.get("locus_tag", ["-"])[0]
        start = int(cds.location.start)
        end = int(cds.location.end)
        strand = cds.location.strand

        categoria = "outros"
        for cat, keywords in categorias.items():
            if any(k in produto for k in keywords):
                categoria = cat
                break

        contagem_funcional[categoria] += 1

        genes_info.append({
            "gene": gene,
            "locus": locus,
            "product": produto,
            "start": start,
            "end": end,
            "strand": strand,
            "category": categoria
        })

    os.makedirs(output_dir, exist_ok=True)

    summary_file = os.path.join(output_dir, "genome_summary.txt")
    with open(summary_file, "w") as f:
        f.write("ANÁLISE GLOBAL DO GENOMA DO BACTERIÓFAGO\n")
        f.write(f"Fago analisado: {phage_name}\n")
        f.write(f"Accession (RefSeq): {accession}\n")
        f.write(f"Tamanho do genoma: {genome_size} bp\n")
        f.write(f"Número total de CDS: {total_cds}\n\n")

        f.write("Distribuição funcional dos genes:\n")
        for cat, count in contagem_funcional.items():
            f.write(f" - {cat}: {count}\n")

    table_file = os.path.join(output_dir, "genome_genes.tsv")
    with open(table_file, "w") as f:
        f.write("gene\tlocus_tag\tproduct\tstart\tend\tstrand\tcategory\n")
        for g in genes_info:
            f.write(
                f"{g['gene']}\t{g['locus']}\t{g['product']}\t"
                f"{g['start']}\t{g['end']}\t{g['strand']}\t{g['category']}\n"
            )

    print("\n--- Análise do genoma concluída ---")
    print(f"Fago: {phage_name}")
    print(f"Accession: {accession}")
    print(f"Tamanho do genoma: {genome_size} bp")
    print(f"Número de CDS: {total_cds}")
    print(f"Resumo salvo em: {summary_file}")
    print(f"Tabela de genes salva em: {table_file}")

    return {
        "phage": phage_name,
        "accession": accession,
        "genome_size": genome_size,
        "total_cds": total_cds,
        "functional_distribution": contagem_funcional
    }
