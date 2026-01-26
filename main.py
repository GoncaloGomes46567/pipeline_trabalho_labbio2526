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
    parser.add_argument("--muscle_path", default="muscle", help="Caminho para o executável muscle.exe")

    args = parser.parse_args()

    
    if not args.email or not args.gene_id:
        args.email, args.gene_id = inputs()

    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)

    gene_info = ncbi.obter_info_gene_ncbi(args.gene_id)
    if not gene_info:
        sys.exit("Erro ao obter informações do gene.")

    gb_file = ncbi.buscar_gene_db(gene_info, args.output)
    if not gb_file:
        sys.exit("Erro no download do GenBank.")

    faa_file = ncbi.extrair_cds_proteina(gb_file, gene_info, args.output)
    if not faa_file:
        sys.exit("Erro na extração da proteína.")


    analisa.analisar_propriedades(faa_file)

    print("\n--- A iniciar o BLAST remoto (isto pode demorar) ---")
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
