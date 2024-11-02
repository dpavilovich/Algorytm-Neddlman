import sys

def read_fasta(file_path):
    
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            sequences.append(line.strip())
    if len(sequences) < 2:
        raise ValueError("Plik FASTA musi zawierać dwie sekwencje.")
    return sequences[0], sequences[1]

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
   
    
    n, m = len(seq1), len(seq2)
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

   
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    
    align1, align2 = '', ''
    i, j = n, m
    while i > 0 and j > 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i - 1][j - 1]
        up_score = score_matrix[i - 1][j]
        left_score = score_matrix[i][j - 1]

        if current_score == diagonal_score + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif current_score == up_score + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    while i > 0:
        align1 = seq1[i - 1] + align1
        align2 = '-' + align2
        i -= 1
    while j > 0:
        align1 = '-' + align1
        align2 = seq2[j - 1] + align2
        j -= 1

    return align1, align2

def calculate_identity_percentage(align1, align2):
    
    matches = sum(1 for a, b in zip(align1, align2) if a == b)
    return (matches / len(align1)) * 100

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Użycie: python nw.py <ścieżka_do_pliku_fasta>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    seq1, seq2 = read_fasta(fasta_path)
    alignment = needleman_wunsch(seq1, seq2)

   
    identity_percentage = calculate_identity_percentage(alignment[0], alignment[1])

    
    with open("alignment_output.txt", 'w') as output_file:
        output_file.write(f"Dopasowanie:\n{alignment[0]}\n{alignment[1]}\n")
        output_file.write(f"Procent identyczności: {identity_percentage:.2f}%\n")
    print("Dopasowanie zapisano w pliku alignment_output.txt")
    print(f"Procent identyczności: {identity_percentage:.2f}%")
