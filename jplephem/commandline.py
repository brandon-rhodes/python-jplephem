import argparse
from .spk import SPK

def main():
    parser = argparse.ArgumentParser(
        prog='python -m jplephem',
        description='Describe an SPK kernel',
        )
    parser.add_argument('path', help='Path to a .bsp SPICE kernel file')
    args = parser.parse_args()
    with open(args.path) as f:
        print(SPK(f))
