# Function to read input file
def read_fasta(file_path):
    """Reads a FASTA file and returns a dictionary with headers as keys and sequences as values."""
    sequences = {} # Initialize an empty dictionary to store sequences
    header = None # Initialize variables to store header and sequence parts
    sequence_parts = []
    with open(file_path, 'r') as file: # Open the FASTA file for reading
        for line in file:  # Iterate over each line in the file
            line = line.strip() # Remove leading and trailing whitespace
            if line.startswith('>'): # Check if line starts with '>'
                if header: # If a header is already present, store the previous sequence
                    sequences[header] = ''.join(sequence_parts).replace(' ', '').upper()
                header = line[1:]    # Extract header (remove '>')
                sequence_parts = [] # Reset sequence parts

            else:
                sequence_parts.append(line)# If the line does not start with '>', append to sequence parts
        if header: # After reaching the end of the file, store the last sequence
            sequences[header] = ''.join(sequence_parts).replace(' ', '').upper()
    return sequences # Return the dictionary of sequences



# Function to find ORFs
def find_orfs(seq, min_orf_length):
    """
    Finds all Open Reading Frames (ORFs) in a DNA sequence for all six reading frames.
    
    Args:
        seq (str): The DNA sequence to search for ORFs.
        min_orf_length (int): The minimum length for an ORF to be considered valid.
    
    Returns:
        list of tuples: Contains information about each ORF found, including frame, start position, and sequence.
    """
    def get_orfs_for_frame(sequence, frame):
        """
        Identifies ORFs in a single reading frame.
        
        Args:
            sequence (str): The DNA sequence to search.
            frame (int): The frame number (0, 1, or 2 for forward; 3, 4, or 5 for reverse complement).
        
        Returns:
            list of tuples: Each tuple contains the start position and the sequence of the ORF.
        """
        start_positions = []  # Track the positions of 'ATG' start codons
        orfs = []  # List to store the ORFs
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == "ATG":  # Start codon
                start_positions.append(i)
            elif codon in ["TAA", "TAG", "TGA"]:  # Stop codons
                # Check each recorded start position to see if it forms a valid ORF with the stop codon
                for start in start_positions:
                    orf_length = i + 3 - start
                    if orf_length >= min_orf_length:
                        orfs.append((start, sequence[start:i+3]))
                start_positions = []  # Reset start positions after finding stop codon
        return orfs
    
    results = []
    # Forward strands: checking frames 0, 1, 2
    for frame in range(3):
        orfs = get_orfs_for_frame(seq, frame)
        for start, orf in orfs:
            results.append((frame + 1, start + 1, orf))  # Store ORF data

    # Reverse complement: generating for frames 3, 4, 5
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_seq = ''.join(complement[nuc] for nuc in seq[::-1])
    for frame in range(3):
        orfs = get_orfs_for_frame(reversed_seq, frame)
        for start, orf in orfs:
            results.append((frame + 4, -start, orf))  # Store reverse ORF data with negative position

    return results

#Function to format ORFs
def format_orf(header, frame, pos, orf):
    """Formats a single ORF in FASTA format with the specified header information."""
    formatted_orf = []
    
    # Split the ORF sequence into codons (triplets)
    for i in range(0, len(orf), 3):
        formatted_orf.append(orf[i:i+3])
    formatted_lines = [' '.join(formatted_orf[i:i+15]) for i in range(0, len(formatted_orf), 15)]

    # Join the formatted lines into a single string with newline characters
    formatted_orf_text = "\n".join(formatted_lines)

    # Create the header string including information like header, frame, position, and length
    orf_header = f">{header} | FRAME = {frame} POS = {pos} LEN = {len(orf)}"

    return f"{orf_header}\n{formatted_orf_text}\n"

def main():
    import os
    file_name = input("Enter the filename of the FASTA file: ")  
    current_directory = os.path.dirname(os.path.abspath(__file__))  
    file_path = os.path.join(current_directory, file_name)  

    min_orf_length = int(input("Enter the minimum ORF length (default is 50): ") or 50)
    
    sequences = read_fasta(file_path)
    for header, seq in sequences.items():
        orfs = find_orfs(seq, min_orf_length)
        for frame, pos, orf in orfs:
            print(format_orf(header, frame, pos, orf))

if __name__ == "__main__":
    main()
