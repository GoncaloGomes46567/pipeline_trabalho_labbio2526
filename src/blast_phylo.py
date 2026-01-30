import os
from collections import defaultdict
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Phylo, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def aplica_blast(fasta_file, output_dir):
    record = SeqIO.read(fasta_file, "fasta")
    print(f"\nA iniciar o BLASTP remoto para {len(record.seq)} aa... (Aguarde)")

    try:
        result_handle = NCBIWWW.qblast(
            "blastp", "nr", record.format("fasta"), hitlist_size=20
        )

        xml_file = os.path.join(output_dir, "blast_results.xml")
        with open(xml_file, "w") as f:
            f.write(result_handle.read())

        result_handle.close()
        return xml_file

    except Exception as e:
        print(f"Erro no BLAST: {e}")
        return None


def processar_blast_e_salvar_hits(xml_file, output_dir, top_hits=10):
    try:
        with open(xml_file) as handle:
            blast_record = NCBIXML.read(handle)

        hits_records = []
        label_counter = defaultdict(int)

        for i, alignment in enumerate(blast_record.alignments[:top_hits]):
            hsp = alignment.hsps[0]
            descricao = alignment.hit_def
            accession = alignment.accession

            # --- proteína ---
            protein = descricao.split("[")[0].lower()
            protein = protein.replace("protein", "").replace("putative", "").strip()
            protein = protein.replace(" ", "_")[:20]

            # --- organismo ---
            organism = "unknown"
            if "[" in descricao and "]" in descricao:
                organism = descricao.split("[")[-1].split("]")[0]
                organism = organism.replace(" ", "_")[:20]

            label = f"{protein}|{organism}"

            label_counter[label] += 1
            if label_counter[label] > 1:
                label = f"{label}_{label_counter[label]}"

            record = SeqRecord(
                Seq(hsp.sbjct),
                id=label,
                description=accession
            )

            hits_records.append(record)

        hits_file = os.path.join(output_dir, "blast_hits.fasta")
        SeqIO.write(hits_records, hits_file, "fasta")

        return hits_file

    except Exception as e:
        print(f"Erro ao processar BLAST: {e}")
        return None


def gerar_arvore_simples(aln_file, output_dir):
    try:
        try:
            alignment = AlignIO.read(aln_file, "fasta")
        except:
            alignment = AlignIO.read(aln_file, "clustal")

        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(alignment)

        tree = DistanceTreeConstructor().nj(dm)

        Phylo.draw_ascii(tree)

        tree_path = os.path.join(output_dir, "arvore.nwk")
        Phylo.write(tree, tree_path, "newick")

        return tree

    except Exception as e:
        print(f"Erro ao gerar a árvore: {e}")
        return None
