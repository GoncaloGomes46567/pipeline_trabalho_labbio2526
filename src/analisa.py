from Bio import SeqIO
from Bio import Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def analisar_propriedades(fasta_file):
    try:
	record = SeqIO.read(fasta_file, "fasta")
	seq_str = str(record.seq)

	analise = ProteinAnalysis(seq_str)
	mw = analise.molecular.weight()
	pi = analise.isoelectric.point()
	instability = analise.instability_index()
	gravy = analise.gravy()

	print(f"\n--- Propriedades Físico-Químicas ({record.id}) ---")
	print(f" - Peso Molecular: {mw: .2f} Da")
	print(f" - Ponto Isoelétrico (pI): {pi:2.f}")
	print(f" - GRAVY: {gravy:.2f}")

	print(" -> Proteína Hidrofóbica" if gravy > 0 else " -> Proteína Hidrofílica")
	return seq_str
    except Exception as e:
	print(f"Erro ao analisar as propriedades: {e}")
	return None

def pesquisar_dominios_cdd(fasta_file):
	pass

