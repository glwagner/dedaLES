if __name__ == '__main__':

from docopt import docopt
from dedalus.tools import logging
from dedalus.tools import post

args = docopt(__doc__)
post.merge_analysis(args['<base_path>'], cleanup=args['--cleanup'])
