import os
import subprocess

def realizar_alinhamento_msa(input_fasta, output_dir, muscle_exe="muscle.exe"):
    output_aln = os.path.join(output_dir, "alinhamento.aln")
    
    print(f"\n--- Iniciando o Alinhamento Múltiplo (MUSCLE v5+) ---")
    
    comando = [
        muscle_exe,
        "-align", input_fasta,   
        "-output", output_aln    
    ]

    try:
        resultado = subprocess.run(comando, capture_output=True, text=True, shell=True)
        
        if resultado.returncode == 0:
            print(f"Alinhamento concluído com sucesso.")
            return output_aln
        else:
            print(f"Erro no MUSCLE: {resultado.stderr}")
            return None
    except Exception as e:
        print(f"Erro ao executar MUSCLE: {e}")
        return None