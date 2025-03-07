from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import *
#
optional_args_dict = {
    **L2D_IsingDynamicsSerializer_optional_args_dict
}
action_args_dict = {
    **L2D_IsingDynamicsSerializer_action_args_dict
}
#
parser = argparse.ArgumentParser(
    description=L2D_IsingDynamicsSerializer_description, 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False
)
#
for k,v in {**optional_args_dict, **action_args_dict}.items():
    parser.add_argument(*k, **v)