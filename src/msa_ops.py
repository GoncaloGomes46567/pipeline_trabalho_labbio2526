import os
import subprocess

import os
import subprocess

def realizar_alinhamento_msa(input_fasta, output_dir, muscle_exe="muscle"):
    """
    Realiza alinhamento múltiplo de sequências usando o MUSCLE.
    """
    output_aln = os.path.join(output_dir, "alinhamento.aln")
    
    print(f"\n--- Iniciando o Alinhamento Múltiplo (MUSCLE) ---")
    
    if not os.path.exists(input_fasta):
        print(f"Erro: Arquivo {input_fasta} não encontrado.")
        return None

    comando = [
        muscle_exe,
        "-in", input_fasta,
        "-out", output_aln
    ]

    try:
        
        resultado = subprocess.run(comando, capture_output=True, text=True)
        
        if resultado.returncode == 0:
            print(f"Alinhamento concluído com sucesso.")
            print(f"Arquivo salvo em: {output_aln}")
            return output_aln
        else:
            print(f"Erro no MUSCLE: {resultado.stderr}")
            return None
            
    except FileNotFoundError:
        print(f"Erro: O executável '{muscle_exe}' não foi encontrado na pasta ou no PATH.")
        return None
    except Exception as e:
        print(f"Erro inesperado: {e}")
        return None