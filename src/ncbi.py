import os
from Bio import Entrez, SeqIO

def setup_entrez(email):
    Entrez.email = email

def procurar_gene(gene_id):
	print(f"A procurar pelo gene {gene_id} na base de dados.")
	try:
		handle = Entrez.esearch(db = "gene", term = gene_id, retmax = 10)
		record = Entrez.read(handle)
		handle.close()
		return record["IdList"]
	except Exception as e:
		print(f"Erro na pesquisa do gene:{e}")
		return []

def buscar_gene_db(gene_id, output_dir):
	try:
		handle_search = Entrez.esearch(db = "nucleotide", term = f"{gene_id}[Gene ID]", retmax = 1)
		record_seach = Entrez.read(handle_search)
		nucleotide_ids = record_seach["IdList"]
		handle_search.close()

		if not nucleotide_ids:
			print(f"Nenhuma sequência encontrada para o gene {gene_id}.")
			return None, None

		seq_id = nucleotide_ids[0]
		print(f"A fazer download da sequência de nucleóitdos (ID : {seq_id})...")

		handle_fetch = Entrez.efetch(db="nucleotide", id = seq_id, rettype = "gb", retmode = "text")
		data = handle_fetch.read()
		handle_fetch.close()

		if not os.path.exists(output_dir):
			os.makedirs(output_dir)

		filename = os.path.join(output_dir, f"gene_{gene_id}.gb")
		with open(filename, "w") as f:
			f.write(data)
		print(f"Arquivo salvo em: {filename}")
		return filename, seq_id

	except Exception as e:
		print(f"Erro ao procurar no Genbank: {e}")
		return None, None

def extrair_cds_proteina(genbank_file, output_dir, gene_target_id="1258817"):
    record = SeqIO.read(genbank_file, "genbank")
    output_faa = os.path.join(output_dir, "protein_sequence.faa")

    print(f"A analisar o registo {record.id} à procura do GeneID: {gene_target_id}...")

    for feature in record.features:
        if feature.type == "CDS":
            db_xrefs = feature.qualifiers.get("db_xref", [])
            is_target = any(f"GeneID:{gene_target_id}" in xref for xref in db_xrefs)

            if is_target:
                translation = None
                if "translation" in feature.qualifiers:
                    translation = feature.qualifiers["translation"][0]
                else:
                    try:
                        translation = str(feature.extract(record.seq).translate(to_stop=True))
                    except Exception as e:
                        print(f"Erro ao traduzir sequência: {e}")
                        continue

                if translation:
                    gene_name = feature.qualifiers.get("gene", ["Desconhecido"])[0]
                    product = feature.qualifiers.get("product", ["Produto desconhecido"])[0]
                    
                    print(f"\n✅ SUCESSO! Encontrado o gene alvo: {gene_name}")
                    print(f"Produto: {product}")
                    
                    with open(output_faa, "w") as f:
                        f.write(f">{gene_name} | GeneID:{gene_target_id} | {product}\n{translation}\n")
                    return output_faa

    print(f"Erro: O GeneID {gene_target_id} não foi encontrado no ficheiro.")
    return None

