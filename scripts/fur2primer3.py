#! /usr/bin/python
import argparse
from pathlib import Path
from Bio import SeqIO # type: ignore


def args_to_dict(params: str)->dict:
    '''The parameter string is split by whitespace and parameters are used as keys and the parametervalues as values for a dictionary. For example primMaxTm=62 would be the key-value pair primMaxTm: 62 in the dictionary.'''
    return {param.split("=")[0]: param.split("=")[1] for param in params.split(" ")}

def remap_keys(params: dict)->dict:
    '''Rename the keys of the dictionary to Primer3 conventions and set defaults.'''
    return {
        "PRIMER_TASK": "generic",
		"PRIMER_PICK_LEFT_PRIMER": "1",
		"PRIMER_PICK_RIGHT_PRIMER": "1",
		"PRIMER_PICK_INTERNAL_OLIGO": params.get("Oligo","1"),
		"PRIMER_MIN_SIZE": params.get("primMinSize", "15"),
		"PRIMER_MAX_SIZE": params.get("primMaxSize","25"),
		"PRIMER_PRODUCT_SIZE_RANGE": f'{params.get("prodMinSize","100")}-{params.get("prodMaxSize", "200")}',
		"PRIMER_MIN_TM": params.get("primMinTm", "58"),
		"PRIMER_OPT_TM": params.get("primOptTm", "60"),
		"PRIMER_MAX_TM": params.get("primMaxTm","62"),
		"PRIMER_INTERNAL_MIN_TM": params.get("inMinTm","63"),
		"PRIMER_INTERNAL_OPT_TM": params.get("inOptTm", "65"),
		"PRIMER_INTERNAL_MAX_TM": params.get("inMaxTm","67")
    }  

def write_result(file:Path, params:dict):
    '''Retrieve the sequences from the fasta, reformat them into Primer3-formatted output using the dictionary.'''
    sequences = SeqIO.parse(file, "fasta")
    outfile=file.with_suffix('.primers.txt')
    with open(outfile, 'w', encoding="utf-8") as f:
        for record in sequences:
            f.writelines([f"{key}={value}\n" for key, value in params.items()])
            f.write(f"SEQUENCE_TEMPLATE={record.seq}\n")
            f.write("=\n")

            

def main():
    '''Run the functions. Takes a fasta and converts it into a Primer3-compatible input file. Do what main does (aka parse arguments, run everything).'''
    #parse the args
    parser = argparse.ArgumentParser(prog="DiPPER2",description='Converting Fastas (including FUR output) to Primer3-compatible output format', epilog="t.wacker2@exeter.ac.uk")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.1.0-beta")
    parser.add_argument('-p', '--parameters', default="primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200", type=str, help="Primer3 primer parameters. Default: primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200. Internal Probe is also picked by default." )
    parser.add_argument('-i', '--input_file', type=Path, required=True, help="The fasta file to be converted. Can be the FUR output file.")

    args = parser.parse_args()
    
    # parse parameters (dictionary)
    params=args_to_dict(args.parameters)

    #remap/ convert parameter keys to Primer3 conventions
    form_param=remap_keys(params)

    #write Primer3 compatible file.
    write_result(args.input_file,form_param)
if __name__ == '__main__':
    main()
 