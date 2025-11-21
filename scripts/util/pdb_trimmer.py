import os
import argparse

def trim_pdb_by_residue_range(input_pdb_path, output_pdb_path, start_res_id, end_res_id, chain_id=None):
    """
    Trims (removes) residues within the specified range, keeping the segments
    that are outside this range. Optionally restricts trimming to a specific chain.
    
    This function implements the inverse of standard trimming: it removes the middle 
    segment (inclusive) and keeps the N-terminal and C-terminal segments.

    Args:
        input_pdb_path (str): Path to the input AlphaFold PDB file.
        output_pdb_path (str): Path where the trimmed PDB will be saved.
        start_res_id (int): The starting residue number (inclusive) of the segment to REMOVE.
        end_res_id (int): The ending residue number (inclusive) of the segment to REMOVE.
        chain_id (str, optional): The ID of the chain to trim (e.g., 'A').
                                  If None, all chains are processed.
    """
    trimmed_lines = []
    
    # PDB format fixed-width column definitions for ATOM records
    # Columns 23-26 contain the Residue Sequence Number (ResID)
    # Column 22 contains the Chain Identifier
    
    try:
        with open(input_pdb_path, 'r') as infile:
            for line in infile:
                # We only care about ATOM, HETATM, and structural records.
                record_type = line[0:6].strip()
                
                if record_type in ['ATOM', 'HETATM']:
                    # Extract the residue number and chain ID
                    res_seq_num_str = line[22:26].strip()
                    current_chain_id = line[21].strip()
                    
                    if not res_seq_num_str:
                        # Skip lines where ResID is missing/unclear
                        continue
                        
                    try:
                        res_seq_num = int(res_seq_num_str)
                    except ValueError:
                        # Handle cases with insertion codes (e.g., 100A).
                        print(f"Warning: Skipping residue with insertion code at line: {line.strip()}")
                        continue
                    
                    # 1. Check if the line belongs to the desired chain (if specified)
                    chain_match = (chain_id is None) or (current_chain_id == chain_id)

                    # 2. Check if the residue is OUTSIDE the trimming range (i.e., we KEEP it)
                    is_outside_range = not (start_res_id <= res_seq_num <= end_res_id)
                    
                    if chain_match and is_outside_range:
                        trimmed_lines.append(line)
                        
                elif record_type == 'TER':
                    # Add TER records, but only if they follow an ATOM/HETATM line that was kept.
                    if trimmed_lines and trimmed_lines[-1].startswith(('ATOM', 'HETATM')):
                        trimmed_lines.append(line)
                        
                elif record_type in ['REMARK', 'MODEL', 'ENDMDL', 'CRYST1', 'MASTER', 'END']:
                    # Keep metadata lines
                    trimmed_lines.append(line)


        # Add the final END record and ensure a clean TER record is present.
        if not trimmed_lines or not trimmed_lines[-1].strip() == 'END':
             # Find the last ATOM line in the kept lines
            last_atom_index = -1
            for i in range(len(trimmed_lines) - 1, -1, -1):
                if trimmed_lines[i].startswith(('ATOM', 'HETATM')):
                    last_atom_index = i
                    break
            
            # If the last kept line is an ATOM, ensure it's followed by a TER record
            if last_atom_index != -1 and not trimmed_lines[last_atom_index].startswith('TER'):
                 last_line = trimmed_lines[last_atom_index]
                 serial = last_line[6:11]
                 resName = last_line[17:20]
                 chainID = last_line[21]
                 resSeq = last_line[22:26]
                 # Create a TER record based on the last atom's numbering
                 ter_line = f"TER   {int(serial):>4}      {resName} {chainID}{resSeq}\n"
                 trimmed_lines.append(ter_line)
                 
            # Ensure the file ends with END
            if not trimmed_lines or not trimmed_lines[-1].strip() == 'END':
                trimmed_lines.append("END\n")

        with open(output_pdb_path, 'w') as outfile:
            outfile.writelines(trimmed_lines)

        print(f"Successfully trimmed PDB file and saved to: {output_pdb_path}")
        print(f"Residues REMOVED in range: {start_res_id} to {end_res_id}")
        print(f"Remaining residues (outside the range) were kept.")
        if chain_id:
            print(f"Chain processed: {chain_id}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_pdb_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Example usage:
    # python trim_pdb.py -i your_alphafold_target.pdb -o trimmed_target.pdb -s 10 -e 150 -c A
    
    parser = argparse.ArgumentParser(description="Trim residues from an AlphaFold PDB file.")
    parser.add_argument('-i', '--input_pdb', required=True, help="Path to the input PDB file.")
    parser.add_argument('-o', '--output_pdb', required=True, help="Path for the output trimmed PDB file.")
    parser.add_argument('-s', '--start_res', type=int, required=True, help="Starting residue number (inclusive) of the segment to REMOVE.")
    parser.add_argument('-e', '--end_res', type=int, required=True, help="Ending residue number (inclusive) of the segment to REMOVE.")
    parser.add_argument('-c', '--chain_id', type=str, default=None, help="Optional: Specify the chain ID to trim (e.g., 'A').")
    
    args = parser.parse_args()

    # The Chain ID often needs to be stripped of whitespace if used in PDB lines
    if args.chain_id:
        args.chain_id = args.chain_id.strip()

    trim_pdb_by_residue_range(
        args.input_pdb,
        args.output_pdb,
        args.start_res,
        args.end_res,
        args.chain_id
    )
