import os
from Bio import Entrez, SeqIO

def setup_entrez:
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
        return[]

def buscar_gene_db(gene_id, output_dir):
    try:
       handle_search: Entrez.esearch(db = "nucleotide", term = f"{gene_id}[Gene ID]", retmax = 1)
       record_seach = Entrez.read(handle_search)
       nucleotide_ids = record_seach["IdList"]
       
       if not nucleotide_ids:
	      print(f"Nenhuma sequência encontrada para o gene {gene_id}.")
	      return None, None
       
       seq_id = nucleotide_ids[0]
       print(f"A fazer download da sequência de nucleóitdos (ID : {seq_id})...")
       
       handle_fetch = Entrez.efetch(db="nucleotide", id = seq_id, rettype = "gb", retmode = "text")
       data = handle_fetch.read()
       handle_fetch.close()

	filename = os.path.join(output_dir, f"gene_{gene_id}.gb")
	with open(filename, "w") as f:
	     f.write(data)
	print(f"Arquivo salvo em: {filename}")
	return filename, seq_id
   except Exception as e:
	print(f"Erro ao procurar no Genbank: {e}")
	return None, None

def extrair_cds_proteina(genbank_file, output_dir, gen_target_name = None):
	record = SeqIO.read(genbank_file, "genbank")
	cds_found = False
	output_faa = os.path.join(output_dir, "protein_sequence.faa")

	for feature in record.features:
		if feature.type == "CDS":
			gene_name = feature.qualifiers.get("gene", ["?"])[0]
			product = feature.qualifiers.get("product",["?"])[0]
			translation = feature.qualifiers.get("traducao",[""])[0]
			match = False
			if gene_target_name:
			   if gene_target_name.lower() in gene_name.lower() or gene_target.lower() in product.lower():
			   match = True
			if match and translation:
			   print(f"\n>>> CDS encontrado com sucesso!")
			   print(f" Gene: {gene_name}")
                           print(f" Produto: {product}")
                           print(f" Tamanho: {len(translation)} aa")
			   with open(output_faa, "w") as f:
			   	.write(f">{gene_name}_{product.replace(" ","-")}\n{translation}\n")
			   print(f" Sequência da proteína salva em: {output_faa})
			   return output_faa
	if not cds_found:
		print("Nenhum CDS com tradução encontrada.")
		return None
