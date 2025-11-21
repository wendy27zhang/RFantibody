import sys
import argparse
import os

def renumber_pdb_residues(input_filename, output_filename, offset):
    """
    Renames the residue sequence numbers (resSeq) in a PDB file by adding a specified offset.
    Only processes ATOM and HETATM records.

    Args:
        input_filename (str): Path to the original PDB file.
        output_filename (str): Path to the output PDB file.
        offset (int): The number to add to the existing resSeq of every residue.
    """
    
    # PDB format uses fixed columns:
    # resSeq is in columns 23-26 (Python slice indices 22:26)
    RESIDUE_NUMBER_START = 22
    RESIDUE_NUMBER_END = 26
    
    # Store the renumbered lines
    modified_lines = []
    
    try:
        with open(input_filename, 'r') as infile:
            print(f"Reading input file: {input_filename}")
            for line in infile:
                # Check if the record is an ATOM or HETATM, which contain coordinates
                record_name = line[0:6].strip()
                
                if record_name in ('ATOM', 'HETATM'):
                    # 1. Extract the current residue number string
                    try:
                        current_resi_str = line[RESIDUE_NUMBER_START:RESIDUE_NUMBER_END].strip()
                        current_resi = int(current_resi_str)
                    except ValueError:
                        # Handle lines where resSeq might be non-numeric or missing
                        modified_lines.append(line)
                        continue
                        
                    # 2. Calculate the new residue number
                    new_resi = current_resi + offset
                    
                    # 3. Format the new number back into a 4-character string
                    # It must be right-justified and padded with spaces to maintain PDB format
                    new_resi_str = f"{new_resi:4d}"
                    
                    if len(new_resi_str) > 4:
                        print(f"Warning: Residue number {new_resi} exceeds 4 digits (max 9999). Residue will be truncated in the output.")

                    # 4. Construct the new line by splicing the string
                    new_line = (
                        line[:RESIDUE_NUMBER_START] + 
                        new_resi_str + 
                        line[RESIDUE_NUMBER_END:]
                    )
                    modified_lines.append(new_line)
                else:
                    # Keep non-coordinate lines (HEADER, REMARK, TER, END, etc.) as is
                    modified_lines.append(line)

        # Write the modified content to the output file
        with open(output_filename, 'w') as outfile:
            outfile.writelines(modified_lines)
            
        print(f"\nSuccessfully renumbered residues.")
        print(f"Offset applied: +{offset}")
        print(f"Output saved to: {output_filename}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        sys.exit(1) # Exit with an error code
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

def main():
    """
    Main function to parse command-line arguments and run the renumbering.
    """
    parser = argparse.ArgumentParser(
        description="Renumber residues in a PDB file by adding a constant offset."
    )
    parser.add_argument(
        "input_pdb",
        help="Path to the input PDB file (e.g., structure.pdb)"
    )
    parser.add_argument(
        "output_pdb",
        help="Path for the renumbered output PDB file (e.g., renumbered.pdb)"
    )
    parser.add_argument(
        "offset",
        type=int,
        help="The integer value to add to every residue number (e.g., 250)"
    )

    args = parser.parse_args()
    
    renumber_pdb_residues(args.input_pdb, args.output_pdb, args.offset)

if __name__ == "__main__":
    main()
