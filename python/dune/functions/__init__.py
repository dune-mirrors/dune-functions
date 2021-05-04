from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dune.common.checkconfiguration import assertHave, ConfigurationError

try:
    assertHave("HAVE_DUNE_FUNCTIONS")
except ConfigurationError:
    raise ImportError("dune module dune-functions was not found")

from .globalbasis import defaultGlobalBasis
from .tree import *

registry = dict()
registry["globalBasis"] = {
        "default" : defaultGlobalBasis
    }

generator = SimpleGenerator("GlobalBasis", "Dune::Python")

def load(includes, typeName, *args):
    includes = includes + ["dune/python/functions/globalbasis.hh"]
    moduleName = "globalbasis_" + hashIt(typeName)
    module = generator.load(includes, typeName, moduleName, *args)
    return module
