import argparse
from src.analise_genomica import analisar_genoma_fago

def main():
    parser = argparse.ArgumentParser(
        description="Análise global do genoma de um bacteriófago"
    )
    parser.add_argument(
        "--genbank",
        required=True,
        help="Ficheiro GenBank do genoma do fago"
    )
    parser.add_argument(
        "--output",
        default="data/genoma_fago",
        help="Pasta de saída"
    )

    args = parser.parse_args()

    analisar_genoma_fago(args.genbank, args.output)

if __name__ == "__main__":
    main()
