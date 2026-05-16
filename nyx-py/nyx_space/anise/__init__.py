from nyx_space import anise as __mod

__all__ = []

for el in dir(__mod):
    if el and el[0] != "_":
        locals()[el] = getattr(__mod, el)
        __all__ += [el]


# from nyx_space.anise import (
#     Aberration,
#     Almanac,
#     MetaAlmanac,
#     MetaFile,
#     LocationDhallSet,
#     LocationDhallSetEntry,
#     LocationDataSet,
#     time,
#     analysis,
#     astro,
#     constants,
#     rotation,
#     utils,
# )

# # __all__ = [
# #     # modules
# #     "analysis",
# #     "astro",
# #     "constants",
# #     "time",
# #     "rotation",
# #     "utils",
# #     # root
# #     "Aberration",
# #     "Almanac",
# #     "MetaAlmanac",
# #     "MetaFile",
# #     "LocationDhallSet",
# #     "LocationDhallSetEntry",
# #     "LocationDataSet",
# #     # functions
# #     "exec_gui",
# #     "__version__",
# #     "__doc__",
# #     "__author__",
# # ]
