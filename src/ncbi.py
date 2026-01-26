import os
from Bio import Entrez, SeqIO


def setup_entrez(email):
    Entrez.email = email

def obter_info_gene_ncbi(gene_id):
    try:
        handle = Entrez.efetch(
            db="gene",
            id=gene_id,
            rettype="xml",
            retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()

        if not record:
            return None

        gene = record[0]
        gene_ref = gene["Entrezgene_gene"]["Gene-ref"]

        symbol = gene_ref.get("Gene-ref_locus", "")
        desc = gene_ref.get("Gene-ref_desc", "")

        organism = gene["Entrezgene_source"]["BioSource"]["BioSource_org"]["Org-ref"]["Org-ref_taxname"]

        return {
            "symbol": symbol,
            "description": desc,
            "organism": organism,
            "gene_id": gene_id
        }

    except Exception as e:
        print(f"Erro ao obter info do gene: {e}")
        return None


def buscar_gene_db(gene_info, output_dir):
    try:
        term = f"{gene_info['symbol']}[Gene] AND {gene_info['organism']}[Organism]"

        handle = Entrez.esearch(
            db="nucleotide",
            term=term,
            retmax=10
        )
        search = Entrez.read(handle)
        handle.close()

        if not search["IdList"]:
            print("Nenhuma sequência nucleotídica encontrada.")
            return None

        seq_id = search["IdList"][0]
        print(f"GenBank encontrado: {seq_id}")

        handle = Entrez.efetch(
            db="nucleotide",
            id=seq_id,
            rettype="gb",
            retmode="text"
        )
        data = handle.read()
        handle.close()

        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"{gene_info['symbol']}.gb")

        with open(filename, "w") as f:
            f.write(data)

        return filename

    except Exception as e:
        print(f"Erro na busca GenBank: {e}")
        return None


def cds_corresponde(feature, gene_info):
    gene_names = feature.qualifiers.get("gene", [])
    locus_tags = feature.qualifiers.get("locus_tag", [])
    products = feature.qualifiers.get("product", [])

    symbol = gene_info["symbol"].lower()

    if symbol in [g.lower() for g in gene_names]:
        return True

    if symbol in [l.lower() for l in locus_tags]:
        return True

    for p in products:
        if symbol in p.lower():
            return True

    return False

def extrair_cds_proteina(genbank_file, gene_info, output_dir):
    record = SeqIO.read(genbank_file, "genbank")
    output_faa = os.path.join(output_dir, "protein_sequence.faa")

    print(f"A analisar {record.id}")
    print(f"Gene alvo: {gene_info['symbol']} ({gene_info['gene_id']})")

    for feature in record.features:
        if feature.type != "CDS":
            continue

        if not cds_corresponde(feature, gene_info):
            continue

        gene_name = feature.qualifiers.get("gene", [gene_info["symbol"]])[0]
        product = feature.qualifiers.get("product", ["Produto desconhecido"])[0]

        if "translation" in feature.qualifiers:
            protein = feature.qualifiers["translation"][0]
        else:
            try:
                protein = str(feature.extract(record.seq).translate(to_stop=True))
            except Exception:
                continue

        print("\nGENE ENCONTRADO!")
        print(f"Gene: {gene_name}")
        print(f"Produto: {product}")

        with open(output_faa, "w") as f:
            f.write(f">{gene_name} | GeneID:{gene_info['gene_id']} | {product}\n")
            f.write(protein + "\n")

        return output_faa

    print("CDS correspondente não encontrada.")
    return None

def obter_proteina_por_geneid(gene_id, output_dir):
    try:
        handle = Entrez.elink(
            dbfrom="gene",
            db="protein",
            id=gene_id
        )
        links = Entrez.read(handle)
        handle.close()

        linksets = links[0].get("LinkSetDb", [])
        if not linksets:
            print("Nenhuma proteína associada ao GeneID.")
            return None

        protein_ids = linksets[0]["Link"]
        protein_id = protein_ids[0]["Id"]

        print(f"Proteína associada encontrada: {protein_id}")

        # Fetch da proteína
        handle = Entrez.efetch(
            db="protein",
            id=protein_id,
            rettype="fasta",
            retmode="text"
        )
        fasta = handle.read()
        handle.close()

        os.makedirs(output_dir, exist_ok=True)
        out = os.path.join(output_dir, "protein_sequence.faa")

        with open(out, "w") as f:
            f.write(fasta)

        return out

    except Exception as e:
        print(f"Erro ao obter proteína: {e}")
        return None
