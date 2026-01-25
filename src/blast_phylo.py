import os
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt


def aplica_blast(fasta_file, output_dir):
	record = SeqIO.read(fasta_file, "fasta")
	print(f"\nA iniciar o BLASTP remoto para {len(record.seq)} aa... (Aguarde)")
	
	try:
		result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"), hitlist_size=20)
		xml_file = os.path.join(output_dir, "blast_results.xml")
		
		with open(xml_file, "w") as f:
			f.write(result_handle.read())
		
		result_handle.close()
		print(f"Resultados do BLAST salvos em: {xml_file}")
		return xml_file
	except Exception as e:
		print(f"Erro no BLAST: {e}")
		return None

def processar_blast_e_salvar_hits(xml_file, output_dir):
	try:
		blast_record = NCBIXML.read(open(xml_file))
		if len(blast_record.alignments) == 0:
			print("Nenhum alinhamento encontrado.")
			return None
		print(f" Foram encontrados {len(blast_record.alignments)} alinhamentos.")
		print(f" Top hit: {blast_record.alignments[0].title[:60]}...")
		hits_fasta = os.path.join(output_dir, "blast_hits.fasta")
		with open(hits_fasta, "w") as f:
			for alignment in blast_record.alignments:
				hsp = alignment.hsps[0]
				clean_id = alignment.accession.replace(" ", "_")
				f.write(f">{clean_id}\n{hsp.sbjct}\n")

		print(f"Sequências dos hits salvas em: {hits_fasta}")
		return hits_fasta
	except Exception as e:
		print(f" Erro ao processar Blast: {e}")
		return None

def gerar_arvore_simples(hits_fasta, output_dir):
	print(f"\n--- Gerando Árvore Filogenética ---")
	try:
		alignment = AlignIO.read(aln_file, "clustal")
		calculator = DistanceCalculator("identity")
		dm = calculator.get_distance(alignment)
		constructor = DistanceTreeConstructor()
		tree = constructor.nj(dm)
		tree_file = os.path.join(output_dir, "arvore.xml")
		Phylo.write(tree, tree_file, "phyloxml")
		plt.figure(figsize =(10,8))
		axes = plt.axes(frameon = False)
		Phylo.draw(tree, axes = axes, do_show = False)
		png_file = os.path.join(output_dir, "arvore_filogenetica.png")
		plt.savefig(png_file)
		print(f"Árvore salva (imagem): {png_file}")
		print(f"Árvore salva(dados): {tree_file}")
	except Exception as e:
		print(f"Erro ao gerar a árvore: {e}")