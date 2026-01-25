import argparse
import os
import sys
from src import ncbi, analisa, blast_phylo, msa_ops


def main():
	parser = argparse.ArgumentParser(description = "Pipeline Genética Completa)")
	parser.add_argument("--email", required = True, help="Email para o NCBI")
	parser.add_argument("--gene_id", required = True, help = "ID do gene")
	parser.add_argument("--xml_file", required=True, help= "Fucgeuri XML do BLAST")
	parser.add_argument("--faa_file", required=True, help = "Ficheiro FASTA (.faa)")
	parser.add_argument("--output", default = "data", help = "Pasta de saída")
	parser.add_argument("--filter_name", default = None, help = "Filtro de nome da proteína")
	parser.add_argument("--clustal_path", default = "clustalw2", help = "Path para o executável do ClustalW. Se instalado no sistema, deixa o padrão.")

	args = parser.parse_args()
	if not os.path.exists(args.output):
		os.makedirs(args.output, exist_ok = True)

	ncbi_ops.setup_entrez(args.email)

	gb_file,  _ = ncbi.buscar_gene_gb(args.gene_id, args.output)
	if not gb_file: sys.exit("erro no download.")
	faa_file = ncbi.extrair_cds_proteina(gb_file, args.output, args.filter_name)
	if not faa_file: sys.exit("Erro na extração")

	analisa.analisar_propriedades(faa_file)

	if xml_file:
		hits_file = blast_phylo.processar_blast_e_salvar_hits(
        xml_file,
        args.output,
        faa_file
    )

	if hits_file:
		aln_file = msa_ops.realizar_alinhamento_clustalw(
            hits_file,
            args.output,
            args.clustal_path
        )

	if aln_file:
            blast_phylo.gerar_arvore_phylo(aln_file, args.output)

if __name__ == "__main__":
	main()
