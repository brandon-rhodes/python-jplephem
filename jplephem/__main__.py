from sys import argv, exit
from .commandline import main
exit(main(argv[1:]) or 0)
