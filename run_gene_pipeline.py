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
    parser.add_argument("--muscle_path", default="muscle", help="Caminho para o executável muscle")

    args = parser.parse_args()

    if not args.email or not args.gene_id:
        args.email, args.gene_id = inputs()

    output_dir = input("Digite a pasta de saída: ")

    ncbi.setup_entrez(args.email)

    os.makedirs(output_dir, exist_ok=True)

    gene_info = ncbi.obter_info_gene_ncbi(args.gene_id)
    if not gene_info:
        sys.exit("Erro ao obter informações do gene.")

    faa_file = ncbi.obter_proteina_por_geneid(args.gene_id, output_dir)
    if not faa_file:
        sys.exit("Erro ao obter proteína a partir do GeneID.")

    analisa.analisar_propriedades(faa_file)

    print("\n--- A iniciar o BLAST remoto (isto pode demorar) ---")
    xml_file = blast_phylo.aplica_blast(faa_file, output_dir)

    if xml_file:
        hits_file = blast_phylo.processar_blast_e_salvar_hits(xml_file, output_dir)

        if hits_file:
            aln_file = msa_ops.realizar_alinhamento_msa(
                hits_file,
                output_dir,
                args.muscle_path
            )

            if aln_file:
                blast_phylo.gerar_arvore_simples(aln_file, output_dir)


if __name__ == "__main__":
    main()
