import sys
from .commandline import main
sys.stdout.write(main(sys.argv[1:]))
sys.exit(0)
