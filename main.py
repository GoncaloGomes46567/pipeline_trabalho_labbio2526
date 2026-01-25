import argparse
import os
import sys
from src import ncbi, analisa, blast_phylo, msa_ops

def inputs():
	email = input("Digite seu email para o NCBI: ")
	gene_id = input("Digite o ID do gene: ")
	return email, gene_id

def main():
    parser = argparse.ArgumentParser(description="Pipeline Genética Completa")
    
    parser.add_argument("--email", help="Email para o NCBI")
    parser.add_argument("--gene_id", help="ID do gene")
    parser.add_argument("--output", default="data", help="Pasta de saída")
    parser.add_argument("--filter_name", default=None, help="Filtro de nome da proteína")
    parser.add_argument("--clustal_path", default="clustalw2", help="Caminho do ClustalW")

    args = parser.parse_args()

    
    if not args.email or not args.gene_id:
        args.email, args.gene_id = inputs()

    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)

    ncbi.setup_entrez(args.email)

    gb_file, _ = ncbi.buscar_gene_db(args.gene_id, args.output)
    if not gb_file:
        sys.exit("erro no download.")
    faa_file = ncbi.extrair_cds_proteina(gb_file, args.output, args.filter_name)
    if not faa_file:
        sys.exit("Erro na extração")

    analisa.analisar_propriedades(faa_file)

    print("\n--- A iniciar o BLAST remoto... Isto pode demorar alguns minutos. ---")
    xml_file = blast_phylo.aplica_blast(faa_file, args.output)
    
    if xml_file:
        hits_file = blast_phylo.processar_blast_e_salvar_hits(xml_file, args.output)

        if hits_file:
         
            aln_file = msa_ops.realizar_alinhamento_msa(
                hits_file, 
                args.output, 
                args.muscle_path
            )

            if aln_file:
                blast_phylo.gerar_arvore_simples(aln_file, args.output)
if __name__ == "__main__":
	main()
