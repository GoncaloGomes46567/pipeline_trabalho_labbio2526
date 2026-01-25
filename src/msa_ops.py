import os
import subprocess

import os
import subprocess

def realizar_alinhamento_clustalw(input_fasta, output_dir, clustalw_exe="clustalw2"):
    output_aln = os.path.join(output_dir, "alinhamento.aln")

    print(f"\n--- Iniciando o Alinhamento Múltiplo (ClustalW) ---")
    
    if not os.path.exists(input_fasta):
        print(f"Erro: Arquivo {input_fasta} não encontrado.")
        return None

    comando = [
        clustalw_exe,
        f"-INFILE={input_fasta}",
        f"-OUTFILE={output_aln}",
        "-OUTPUT=CLUSTAL"
    ]

    try:
        resultado = subprocess.run(comando, capture_output=True, text=True)
        
        if resultado.returncode == 0:
            print(f"Alinhamento concluído com sucesso.")
            print(f"Arquivo salvo em: {output_aln}")
            return output_aln
        else:
            print(f"Erro no ClustalW: {resultado.stderr}")
            return None
            
    except FileNotFoundError:
        print(f"Erro: O executável '{clustalw_exe}' não foi encontrado.")
        print("Certifique-se de que o ClustalW2 está instalado e no PATH do Windows.")
        return None
    except Exception as e:
        print(f"Erro inesperado: {e}")
        return None