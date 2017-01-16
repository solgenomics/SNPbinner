'''Sets up the main argparser and imports and sets up scripts in the program.'''

import __init__ as snpbinner
import sys
import argparse

def main():
    program_dict = {}
    for module in snpbinner.__all__:
        program_dict[module] = getattr(snpbinner,module)

    main_parser = argparse.ArgumentParser(description="The snpbinner package can be run in three ways. It can be run as an executable directly from the commandline, run as a python program with ``$python snpbinner`` or it can be imported as a python module :class:`snpbinner` and used with other python scripts.")
    prog_sub = main_parser.add_subparsers(title="Program",dest="program",help="Specifies which program in the snpbinner package should be run.")

    program_run_dict = {name:program_dict[name].run for name in program_dict if hasattr(program_dict[name],"run")}
    program_parser_dict = {name:program_dict[name].parser(prog_sub.add_parser,name) for name in program_dict if name in program_run_dict}

    args = main_parser.parse_args(sys.argv[1:])

    program_to_run = args.program
    del args.program
    arg_dict = vars(args)
    arg_dict_to_pass = {key:arg_dict[key] for key in arg_dict if key!="program" and arg_dict[key]!=None}
    program_run_dict[program_to_run](**arg_dict_to_pass)

if __name__ == '__main__':
    main()