import argparse
from .daf import DAF
from .spk import SPK

def main(args):
    parser = argparse.ArgumentParser(
        prog='python -m jplephem',
        description='Describe an SPK kernel',
        )
    parser.add_argument('path', help='Path to a .bsp SPICE kernel file')
    args = parser.parse_args(args)
    with open(args.path, 'rb') as f:
        return(str(SPK(DAF(f))))
